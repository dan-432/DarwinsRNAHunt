#!/usr/bin/env python3
"""
Grid-test mmseqs2 clustering parameters for cluster size and purity.

Goal: find parameter combinations that produce clusters which are large
(more sequence variance) while still being mostly composed of sequences
with real structural signal against a reference CM — wide enough to be
useful, not so wide that alignment/CM-building quality suffers.

For each parameter combination, mmseqs clustering is followed by extraction
of one FASTA per cluster. Each cluster FASTA is either searched as-is or
trimmed to its conserved window first, and cmsearch is then run on that
cluster FASTA.

For each parameter combination, reports (over clusters at or above
--min-cluster-size):
  n_clusters_total       clusters produced before size filtering
  n_clusters_passing     clusters at or above --min-cluster-size
  mean_cluster_size      mean size of passing clusters
  mean_purity            mean, across passing clusters, of the fraction
                          of members with a significant cmsearch hit
  weighted_purity         significant members / total members, pooled
                          across all passing clusters (sensitive to a
                          few large low-purity clusters, unlike mean_purity)

Usage:
    conda run -n domain-taxonomy python test_mmseqs_params.py \
        --input flanks.fasta \
        --cm RF00050.cm \
        --outdir results/ \
        --min-cluster-size 5

Requires mmseqs2 and Infernal (cmsearch) on PATH. Standard library only.
"""

import argparse
import csv
import itertools
import re
import subprocess
from collections import defaultdict
from pathlib import Path
from statistics import mean

# ---------------------------------------------------------------------------
# Parameter grid — edit to control which combinations get tested.
# Clustering is cheap relative to cmsearch, so this can be widened freely.
# min_seq_id is the main lever for cluster size.
# ---------------------------------------------------------------------------
PARAM_GRID = {
    "sensitivity": [7.5],
    "min_seq_id": [0.2, 0.4, 0.8],
    "coverage_threshold": [0.1, 0.3, 0.5],
    "cov_mode": [1, 3],
    "cluster_mode": [0, 1],
    "trim_to_window": [True]
}


def param_combinations(grid):
    keys = list(grid.keys())
    for values in itertools.product(*(grid[k] for k in keys)):
        yield dict(zip(keys, values))


def run_cmsearch(input_fasta, cm_model, tblout):
    """Run cmsearch on one cluster FASTA and return its tblout path."""
    tblout.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "conda", "run", "-n", "rna-motif", "cmsearch",
        "--tblout", str(tblout),
        "--noali",
        "--cpu", "16",
        str(cm_model), str(input_fasta),
    ]
    subprocess.run(cmd, check=True)
    return tblout


def parse_best_bitscores(tblout, evalue_cutoff):
    """Parse a cmsearch --tblout file into {seq_id: best_bitscore} for hits
    passing the given E-value cutoff. A sequence can appear in multiple
    hits (e.g. multiple domains); only the best-scoring hit is kept."""
    best = {}
    with open(tblout) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = re.split(r"\s+", line.strip())
            seq_id = fields[0]
            score = float(fields[14])
            evalue = float(fields[15])
            if evalue <= evalue_cutoff:
                if seq_id not in best or score > best[seq_id]:
                    best[seq_id] = score
    return best


def run_mmseqs_cluster(input_fasta, workdir, params):
    """Run mmseqs easy-cluster with the given parameters.

    Returns the path to the cluster TSV (representative, member) it
    produces.
    """
    workdir.mkdir(parents=True, exist_ok=True)
    prefix = workdir / "cluster"
    tmp = workdir / "tmp"

    cmd = [
        "conda", "run", "-n", "rna-motif", "mmseqs", "easy-cluster",
        str(input_fasta), str(prefix), str(tmp),
        "--min-seq-id", str(params["min_seq_id"]),
        "-c", str(params["coverage_threshold"]),
        "--cov-mode", str(params["cov_mode"]),
        "--cluster-mode", str(params["cluster_mode"]),
        "--threads", "16",
    ]

    subprocess.run(cmd)

    cluster_tsv = Path(f"{prefix}_cluster.tsv")
    if not cluster_tsv.exists():
        raise FileNotFoundError(f"Expected mmseqs output not found: {cluster_tsv}")
    return cluster_tsv


def run_mmseqs_search(input_fasta, workdir, params):
    """Run an all-vs-all nucleotide search and return its alignment table."""
    db_dir = workdir / "search_db"
    db_dir.mkdir(parents=True, exist_ok=True)
    input_db = db_dir / "input"
    search_db = db_dir / "search"
    tmp = workdir / "search_tmp"
    alignments = workdir / "alignments.m8"

    subprocess.run([
        "conda", "run", "-n", "rna-motif", "mmseqs", "createdb",
        str(input_fasta), str(input_db),
    ])
    subprocess.run([
        "conda", "run", "-n", "rna-motif", "mmseqs", "search",
        str(input_db), str(input_db), str(search_db), str(tmp),
        "--min-seq-id", str(params["min_seq_id"]),
        "-c", str(params["coverage_threshold"]),
        "--cov-mode", str(params["cov_mode"]),
        "--search-type", "3",
        "--alignment-mode", "3",
        "-a",
        "--threads", "16",
    ], check=True)
    subprocess.run([
        "conda", "run", "-n", "rna-motif", "mmseqs", "convertalis",
        str(input_db), str(input_db), str(search_db), str(alignments),
        "--format-output",
        "query,target,qstart,qend,tstart,tend,qlen,tlen,fident,bits,evalue",
    ])
    return alignments


def read_fasta(path):
    """Return {sequence ID: sequence} using the first FASTA header token."""
    sequences = {}
    sequence_id = None
    chunks = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if sequence_id is not None:
                    sequences[sequence_id] = "".join(chunks)
                sequence_id = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if sequence_id is not None:
        sequences[sequence_id] = "".join(chunks)
    return sequences


def write_cluster_fastas(input_fasta, clusters, workdir, min_cluster_size):
    """Write FASTAs only for clusters meeting the minimum size."""
    sequences = read_fasta(input_fasta)
    cluster_fastas = {}
    for index, (representative, members) in enumerate(clusters.items()):
        if len(members) < min_cluster_size:
            continue
        cluster_fasta = workdir / "clusters" / f"cluster_{index:04d}.fasta"
        cluster_fasta.parent.mkdir(parents=True, exist_ok=True)
        with open(cluster_fasta, "w") as fh:
            for member in members:
                if member not in sequences:
                    raise KeyError(f"Sequence {member!r} from {representative!r} is missing")
                fh.write(f">{member}\n{sequences[member]}\n")
        cluster_fastas[representative] = cluster_fasta
    return cluster_fastas


def search_clusters(input_fasta, cm_model, clusters, run_dir, params, evalue,
                    min_cluster_size):
    """Optionally trim and cmsearch clusters meeting the minimum size."""
    cluster_fastas = write_cluster_fastas(
        input_fasta, clusters, run_dir, min_cluster_size
    )
    bitscores = {}
    searchable_clusters = {}
    alignment_path = run_dir / "alignments.m8"
    for representative, cluster_fasta in cluster_fastas.items():
        search_fasta = cluster_fasta
        if params["trim_to_window"]:
            trimmed = cluster_fasta.with_name(f"{cluster_fasta.stem}.trimmed.fasta")
            subprocess.run([
                "python", "workflow/src/trim_to_conserved_window.py",
                "--fasta", str(cluster_fasta),
                "--alignments", str(alignment_path),
                "--output", str(trimmed),
                "--min-length", "30",
                "--max-length", "1000",
                "--padding", "30",
            ], check=True)
            search_fasta = trimmed

            if not read_fasta(search_fasta):
                continue

        tblout = run_cmsearch(
            search_fasta, cm_model,
            run_dir / "cmsearch" / f"{cluster_fasta.stem}.tblout",
        )
        bitscores.update(parse_best_bitscores(tblout, evalue))
        searchable_clusters[representative] = clusters[representative]
    return bitscores, searchable_clusters


def parse_clusters(cluster_tsv):
    """Parse an mmseqs cluster TSV (representative<TAB>member per line)
    into {representative: [member_ids]}."""
    clusters = defaultdict(list)
    with open(cluster_tsv, newline="") as fh:
        for rep, member in csv.reader(fh, delimiter="\t"):
            clusters[rep].append(member)
    return clusters


def score_clusters(clusters, bitscores, min_cluster_size, min_bitscore):
    """Compute size and purity for clusters at or above min_cluster_size.

    A member counts as "significant" if it has a cmsearch hit passing the
    E-value cutoff (i.e. present in `bitscores`) and, if min_bitscore is
    set, its best bitscore also meets that threshold.
    """
    passing = []
    for members in clusters.values():
        if len(members) < min_cluster_size:
            continue
        significant = 0
        for member in members:
            score = bitscores.get(member)
            if score is None:
                continue
            if min_bitscore is not None and score < min_bitscore:
                continue
            significant += 1
        passing.append({
            "size": len(members),
            "significant": significant,
            "purity": significant / len(members),
        })
    return passing


def summarize(clusters, passing):
    if not passing:
        return {
            "n_clusters_total": len(clusters),
            "n_clusters_passing": 0,
            "mean_cluster_size": "",
            "mean_purity": "",
            "weighted_purity": "",
        }

    total_members = sum(c["size"] for c in passing)
    total_significant = sum(c["significant"] for c in passing)
    return {
        "n_clusters_total": len(clusters),
        "n_clusters_passing": len(passing),
        "mean_cluster_size": round(mean(c["size"] for c in passing), 2),
        "mean_purity": round(mean(c["purity"] for c in passing), 3),
        "weighted_purity": round(total_significant / total_members, 3),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path,
                         help="Input FASTA of flanking regions to cluster")
    parser.add_argument("--cm", required=True, type=Path,
                         help="Path to RF00050.cm")
    parser.add_argument("--outdir", required=True, type=Path,
                         help="Directory for per-run outputs and summary.csv")
    parser.add_argument("--evalue", type=float, default=0.01,
                         help="cmsearch E-value cutoff for a significant hit (default 0.01)")
    parser.add_argument("--min-bitscore", type=float, default=None,
                         help="Optional additional bitscore floor for a member to count "
                              "as 'real structural signal' (e.g. the family's GA cutoff). "
                              "Default: rely on --evalue alone.")
    parser.add_argument("--min-cluster-size", type=int, default=3,
                         help="Minimum cluster size to include in summary stats (default 3)")
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    summary_path = args.outdir / "summary.csv"
    fieldnames = list(PARAM_GRID.keys()) + [
        "n_clusters_total", "n_clusters_passing",
        "mean_cluster_size", "mean_purity", "weighted_purity",
    ]

    with open(summary_path, "w", newline="") as summary_fh:
        writer = csv.DictWriter(summary_fh, fieldnames=fieldnames)
        writer.writeheader()

        for i, params in enumerate(param_combinations(PARAM_GRID)):
            run_dir = args.outdir / f"run_{i:03d}"
            print(f"[{i}] {params}")

            try:
                alignment_path = None
                if params["trim_to_window"]:
                    alignment_path = run_mmseqs_search(args.input, run_dir, params)
                cluster_tsv = run_mmseqs_cluster(args.input, run_dir, params)
                clusters = parse_clusters(cluster_tsv)
                bitscores, searchable_clusters = search_clusters(
                    args.input,
                    args.cm,
                    clusters,
                    run_dir,
                    params,
                    args.evalue,
                    args.min_cluster_size,
                )
                passing = score_clusters(
                    searchable_clusters, bitscores,
                    args.min_cluster_size, args.min_bitscore
                )
                row = {**params, **summarize(searchable_clusters, passing)}
            except subprocess.CalledProcessError as e:
                print(f"  FAILED: {e}")
                row = {**params, "n_clusters_total": "", "n_clusters_passing": "",
                       "mean_cluster_size": "", "mean_purity": "", "weighted_purity": ""}

            writer.writerow(row)
            summary_fh.flush()

    print(f"\nDone. Summary written to {summary_path}")


if __name__ == "__main__":
    main()