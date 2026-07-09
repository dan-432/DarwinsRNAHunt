#!/usr/bin/env python3
"""
analyse_cluster_spread.py

Computes per-cluster distribution statistics (sequence counts, genome counts,
prevalence across the available genome set) for a family's flanking-CDS
mmseqs2 clusters, and emits a genome x cluster presence/absence matrix for
downstream phylogenetic analysis (see analyse_cluster_phylogeny.py).
"""

import argparse
import json
import logging
from collections import defaultdict

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--clusters", required=True,
                    help="flanking_cds_clusters.tsv (cols: cluster_id, sequence_id)")
    p.add_argument("--fasta", required=True,
                    help="combined_flanking_cds.faa (headers: genome_accession|protein_id|CDS|...)")
    p.add_argument("--metadata", required=True,
                    help="genomes_available_metadata.json (checkpoint output; key 'accessions')")
    p.add_argument("--output_report", required=True)
    p.add_argument("--output_presence_matrix", required=True)
    p.add_argument("--output_summary", required=True)
    return p.parse_args()


def load_clusters(path):
    """cluster_id -> list of sequence_ids"""
    clusters = defaultdict(list)
    with open(path) as f:
        for line_num, line in enumerate(f, 1):
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                log.warning("Skipping malformed line %d in %s: %r", line_num, path, line)
                continue
            cluster_id, seq_id = parts[0], parts[1]
            clusters[cluster_id].append(seq_id)
    return clusters


def load_seq_to_genome(fasta_path):
    """sequence_id (fasta header token) -> genome_accession (first '|'-delimited field)"""
    seq_to_genome = {}
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                header = line[1:].strip()
                seq_id = header.split()[0]
                genome = seq_id.split("|")[0]
                seq_to_genome[seq_id] = genome
    return seq_to_genome


def load_available_genomes(metadata_path):
    with open(metadata_path) as f:
        data = json.load(f)
    accessions = data.get("accessions", [])
    if not accessions:
        log.warning("No 'accessions' key (or empty list) found in %s", metadata_path)
    return sorted(set(accessions))


def main():
    args = parse_args()

    clusters = load_clusters(args.clusters)
    seq_to_genome = load_seq_to_genome(args.fasta)
    available_genomes = load_available_genomes(args.metadata)
    n_available = len(available_genomes)

    log.info("Loaded %d clusters, %d sequences (fasta), %d available genomes",
              len(clusters), len(seq_to_genome), n_available)

    missing = 0
    report_rows = []
    presence = {}

    for cluster_id, seq_ids in clusters.items():
        genomes_in_cluster = set()
        for seq_id in seq_ids:
            genome = seq_to_genome.get(seq_id)
            if genome is None:
                missing += 1
                continue
            genomes_in_cluster.add(genome)

        n_sequences = len(seq_ids)
        n_genomes = len(genomes_in_cluster)
        prevalence = n_genomes / n_available if n_available else 0.0
        copies_per_genome = n_sequences / n_genomes if n_genomes else 0.0

        presence[cluster_id] = genomes_in_cluster
        report_rows.append({
            "cluster_id": cluster_id,
            "n_sequences": n_sequences,
            "n_genomes": n_genomes,
            "prevalence": round(prevalence, 4),
            "copies_per_genome": round(copies_per_genome, 3),
        })

    if missing:
        log.warning("%d sequence IDs from %s were not found in %s (excluded from genome counts)",
                    missing, args.clusters, args.fasta)

    report_rows.sort(key=lambda r: r["n_genomes"], reverse=True)

    # --- spread report ---
    with open(args.output_report, "w") as f:
        f.write("cluster_id\tn_sequences\tn_genomes\tprevalence\tcopies_per_genome\n")
        for r in report_rows:
            f.write(f"{r['cluster_id']}\t{r['n_sequences']}\t{r['n_genomes']}\t"
                    f"{r['prevalence']}\t{r['copies_per_genome']}\n")

    # --- presence/absence matrix: clusters x ALL available genomes ---
    cluster_ids_sorted = [r["cluster_id"] for r in report_rows]
    with open(args.output_presence_matrix, "w") as f:
        f.write("cluster_id\t" + "\t".join(available_genomes) + "\n")
        for cid in cluster_ids_sorted:
            genomes_present = presence[cid]
            row = [cid] + ["1" if g in genomes_present else "0" for g in available_genomes]
            f.write("\t".join(row) + "\n")

    # --- summary ---
    n_clusters = len(clusters)
    n_singletons = sum(1 for r in report_rows if r["n_genomes"] <= 1)
    top_n = report_rows[:10]

    with open(args.output_summary, "w") as f:
        f.write(f"Total clusters: {n_clusters}\n")
        f.write(f"Total available genomes: {n_available}\n")
        f.write(f"Singleton-genome clusters: {n_singletons}\n")
        if missing:
            f.write(f"WARNING: {missing} sequence IDs missing from fasta headers\n")
        f.write("\nTop 10 clusters by genome spread:\n")
        f.write("cluster_id\tn_genomes\tn_sequences\tprevalence\n")
        for r in top_n:
            f.write(f"{r['cluster_id']}\t{r['n_genomes']}\t{r['n_sequences']}\t{r['prevalence']}\n")

    log.info("Done. Wrote %s, %s, %s",
              args.output_report, args.output_presence_matrix, args.output_summary)


if __name__ == "__main__":
    main()
