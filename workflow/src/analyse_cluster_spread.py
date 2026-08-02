#!/usr/bin/env python3
"""
analyse_cluster_spread.py

Computes per-cluster distribution statistics (sequence counts, genome counts,
prevalence across the available genome set) for a family's flanking-CDS
mmseqs2 clusters, and emits a genome x cluster presence/absence matrix for
downstream phylogenetic analysis (see analyse_cluster_phylogeny.py).

Also folds in per-cluster positional makeup (upstream/downstream/contained
balance, same-strand fraction, fully-intergenic fraction) from the
flanking-CDS position table produced by extract_flanking_cds.py, so you can
see e.g. "cluster_0042 is mostly the immediate 5' same-strand neighbor" at a
glance.
"""

import argparse
import csv
import json
import logging
from collections import defaultdict, Counter
from DarwinsRNAHunt.motif import PositionRecord

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--clusters", required=True,
                    help="flanking_cds_clusters.tsv (cols: cluster_id, sequence_id)")
    p.add_argument("--fasta", required=True,
                    help="combined_flanking_cds.faa (headers: genome_accession|protein_id|Protein|...)")
    p.add_argument("--position_table", required=True,
                    help="combined_flanking_cds_positions.tsv (cols: protein_id, position, rank, "
                         "bp_distance_from_motif, cds_strand, motif_strand, same_strand, seqid, "
                         "cds_start, cds_end, motif_fully_intergenic)")
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


def load_position_table(path):
    """
    Args:
        path to positions table (as written by extract_flanking_cds.py via CDSPosition.to_row)
    Returns:
        dict protein_id (fasta header string, including leading '>') -> CDSPosition
    """
    positions = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            position = PositionRecord.from_row(row)
            positions[row["protein_id"]] = position
    return positions


def load_available_genomes(metadata_path):
    """
    Args:
        metadata_path str: path to file with genome accessions used in analysis
    Returns:
        sorted set of accession ids from file"""
    with open(metadata_path) as f:
        data = json.load(f)
    accessions = data.get("accessions", [])
    if not accessions:
        log.warning("No 'accessions' key (or empty list) found in %s", metadata_path)
    return sorted(set(accessions))


def summarise_cluster_positions(seq_ids, positions):
    """
    Given a cluster's sequence IDs and the protein_id -> PositionRecord lookup,
    return a dict of positional summary stats for that cluster.
    Sequences with no matching position-table entry are excluded from these
    stats (counted separately as n_missing_position).

    Note: median_abs_bp_distance only considers upstream/downstream (disjoint)
    rows. For 'contained' rows, .distance is a coordinate-gap value that can
    be misleading as an overlap magnitude once one interval nests inside the
    other (see gap_to() in sequence.py) -- so it's excluded here rather than
    silently mixed into a "distance" statistic it doesn't represent.
    """
    position_sides = []    # "upstream" / "downstream" / "contained"
    position_labels = []   # "upstream:1", "downstream:2", "contained:0" etc.
    same_strand_flags = []
    intergenic_flags = []
    disjoint_distances = []
    n_missing = 0

    for seq_id in seq_ids:
        pos = positions.get(">" + seq_id)
        if pos is None:
            n_missing += 1
            continue
        position_sides.append(pos.position.value)
        position_labels.append(f"{pos.position.value}:{pos.rank}")
        same_strand_flags.append(pos.same_strand)
        intergenic_flags.append(pos.motif_fully_intergenic)
        if pos.position.value != "contained":
            disjoint_distances.append(pos.distance)

    n_upstream = sum(1 for s in position_sides if s == "upstream")
    n_downstream = sum(1 for s in position_sides if s == "downstream")
    n_contained = sum(1 for s in position_sides if s == "contained")
    n_placed = n_upstream + n_downstream
    frac_upstream = n_upstream / n_placed if n_placed else None
    frac_contained = n_contained / len(position_sides) if position_sides else None

    dominant_position, dominant_count = (None, 0)
    if position_labels:
        dominant_position, dominant_count = Counter(position_labels).most_common(1)[0]
    dominant_position_frac = dominant_count / len(position_labels) if position_labels else None

    frac_same_strand = (sum(same_strand_flags) / len(same_strand_flags)
                        if same_strand_flags else None)
    frac_fully_intergenic = (sum(intergenic_flags) / len(intergenic_flags)
                             if intergenic_flags else None)

    median_abs_distance = None
    if disjoint_distances:
        abs_sorted = sorted(abs(d) for d in disjoint_distances)
        mid = len(abs_sorted) // 2
        median_abs_distance = (abs_sorted[mid] if len(abs_sorted) % 2
                               else (abs_sorted[mid - 1] + abs_sorted[mid]) / 2)

    return {
        "n_missing_position": n_missing,
        "dominant_relative_position": dominant_position if dominant_position else "NA",
        "dominant_relative_position_frac": (round(dominant_position_frac, 3)
                                            if dominant_position_frac is not None else ""),
        "frac_upstream": round(frac_upstream, 3) if frac_upstream is not None else "",
        "frac_contained": round(frac_contained, 3) if frac_contained is not None else "",
        "frac_same_strand": round(frac_same_strand, 3) if frac_same_strand is not None else "",
        "frac_fully_intergenic": round(frac_fully_intergenic, 3) if frac_fully_intergenic is not None else "",
        "median_abs_bp_distance": (round(median_abs_distance, 1)
                                   if median_abs_distance is not None else ""),
    }


def main():
    args = parse_args()

    clusters = load_clusters(args.clusters)
    seq_to_genome = load_seq_to_genome(args.fasta)
    positions = load_position_table(args.position_table)
    available_genomes = load_available_genomes(args.metadata)
    n_available = len(available_genomes)

    log.info("Loaded %d clusters, %d sequences (fasta), %d position records, %d available genomes",
              len(clusters), len(seq_to_genome), len(positions), n_available)

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

        position_summary = summarise_cluster_positions(seq_ids, positions)

        presence[cluster_id] = genomes_in_cluster
        row = {
            "cluster_id": cluster_id,
            "n_sequences": n_sequences,
            "n_genomes": n_genomes,
            "prevalence": round(prevalence, 4),
            "copies_per_genome": round(copies_per_genome, 3),
        }
        row.update(position_summary)
        report_rows.append(row)

    if missing:
        log.warning("%d sequence IDs from %s were not found in %s (excluded from genome counts)",
                    missing, args.clusters, args.fasta)

    n_missing_position_total = sum(r["n_missing_position"] for r in report_rows)
    if n_missing_position_total:
        log.warning("%d sequence IDs across all clusters had no entry in %s "
                    "(excluded from positional stats)",
                    n_missing_position_total, args.position_table)

    report_rows.sort(key=lambda r: r["n_genomes"], reverse=True)

    report_columns = [
        "cluster_id", "n_sequences", "n_genomes", "prevalence", "copies_per_genome",
        "dominant_relative_position", "dominant_relative_position_frac", "frac_upstream",
        "frac_contained", "frac_same_strand", "frac_fully_intergenic", "median_abs_bp_distance",
        "n_missing_position",
    ]

    # --- spread report ---
    with open(args.output_report, "w") as f:
        f.write("\t".join(report_columns) + "\n")
        for r in report_rows:
            f.write("\t".join(str(r[c]) for c in report_columns) + "\n")

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
        if n_missing_position_total:
            f.write(f"WARNING: {n_missing_position_total} sequence IDs missing from position table\n")
        f.write("\nTop 10 clusters by genome spread:\n")
        f.write("cluster_id\tn_genomes\tn_sequences\tprevalence\tdominant_relative_position\t"
                "frac_upstream\tfrac_contained\tfrac_same_strand\tfrac_fully_intergenic\n")
        for r in top_n:
            f.write(f"{r['cluster_id']}\t{r['n_genomes']}\t{r['n_sequences']}\t{r['prevalence']}\t"
                    f"{r['dominant_relative_position']}\t{r['frac_upstream']}\t{r['frac_contained']}\t"
                    f"{r['frac_same_strand']}\t{r['frac_fully_intergenic']}\n")

    log.info("Done. Wrote %s, %s, %s",
              args.output_report, args.output_presence_matrix, args.output_summary)


if __name__ == "__main__":
    main()