"""
analyse_cluster_phylogeny.py

For each cluster in a genome x cluster presence/absence matrix, evaluate how
the cluster's presence is distributed across a genome phylogeny: is it
monophyletic, how "pure" is its MRCA clade, and how many independent
gain/loss events are required to explain the pattern under Fitch (unordered)
and Dollo (single-origin) parsimony. Also emits an iTOL DATASET_BINARY
annotation file for the largest clusters so the tree can be visually
inspected in iTOL.
"""

import argparse
import logging
from itertools import cycle

from ete3 import Tree

from DarwinsRNAHunt.gtdb_access import normalize_accession, get_assembly_to_gtdb_dict, ncbi_accessions_to_gtdb
from DarwinsRNAHunt.gtdb_taxonomy_analysis import trim_tree_to_taxa

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

ITOL_PALETTE = [
    "#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231",
    "#911eb4", "#46f0f0", "#f032e6", "#bcf60c", "#fabebe",
]


def load_presence_matrix(path):
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        genomes = header[1:]
        matrix = {}
        for line in f:
            parts = line.rstrip("\n").split("\t")
            cluster_id, values = parts[0], parts[1:]
            matrix[cluster_id] = dict(zip(genomes, values))
    return genomes, matrix


def build_genome_to_leaf(genomes, assembly_to_gtdb_dict, tree_leaf_names):
    """Resolve each matrix genome accession to exactly one tree leaf name.

    assembly_to_gtdb_dict maps normalized NCBI assembly accessions to a LIST
    of GTDB representative accessions (can be >1 if GTDB reclassified a
    genome across releases). This picks whichever of those candidates is
    actually a leaf in the (already pruned) tree.

    Returns:
        genome_to_leaf: dict genome -> leaf name, or None if unresolved
        unmatched: list of genomes with no candidate present in the tree
    """
    genome_to_leaf = {}
    unmatched = []

    for g in genomes:
        candidates = assembly_to_gtdb_dict.get(normalize_accession(g), [])
        in_tree = [c for c in candidates if c in tree_leaf_names]

        if len(in_tree) == 1:
            genome_to_leaf[g] = in_tree[0]
        elif len(in_tree) > 1:
            log.warning("%s maps to %d tree leaves %s - using first", g, len(in_tree), in_tree)
            genome_to_leaf[g] = in_tree[0]
        else:
            genome_to_leaf[g] = None
            unmatched.append(g)

    if unmatched:
        log.warning("%d/%d matrix genomes had no matching tree leaf (showing up to 10): %s",
                    len(unmatched), len(genomes), unmatched[:10])

    matched_leaf_names = {v for v in genome_to_leaf.values() if v is not None}
    tree_only = tree_leaf_names - matched_leaf_names
    if tree_only:
        log.warning("%d tree leaves have no corresponding entry in the presence matrix "
                    "(showing up to 10): %s", len(tree_only), list(tree_only)[:10])

    return genome_to_leaf, unmatched


def fitch_transitions(tree, presence_leaf_names):
    """Minimum number of state changes for a binary character under Fitch
    parsimony, evaluated over the whole (already-pruned) tree."""
    for leaf in tree.get_leaves():
        leaf.add_feature("_fitch", {1} if leaf.name in presence_leaf_names else {0})

    transitions = 0
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            continue
        child_sets = [c._fitch for c in node.children]
        inter = set.intersection(*child_sets)
        if inter:
            node.add_feature("_fitch", inter)
        else:
            node.add_feature("_fitch", set.union(*child_sets))
            transitions += 1
    return transitions


def dollo_losses(node, presence_leaf_names):
    """Count minimal independent loss events under Dollo parsimony (single
    gain assumed at `node`, the MRCA of all presences)."""
    leaf_names = {l.name for l in node.get_leaves()}
    if leaf_names.isdisjoint(presence_leaf_names):
        return 1  # whole subtree is one loss event; don't recurse further
    if node.is_leaf():
        return 0
    return sum(dollo_losses(child, presence_leaf_names) for child in node.children)


def analyse_cluster(tree, presence_leaf_names):
    if len(presence_leaf_names) < 2:
        return {
            "is_monophyletic": "NA", "clade_type": "NA", "clade_purity": "NA",
            "fitch_transitions": "NA", "dollo_transitions": "NA",
        }

    is_mono, clade_type, _ = tree.check_monophyly(
        values=presence_leaf_names, target_attr="name", ignore_missing=True
    )
    mrca = tree.get_common_ancestor(list(presence_leaf_names))
    clade_purity = round(len(presence_leaf_names) / len(mrca.get_leaves()), 4)

    fitch = fitch_transitions(tree, presence_leaf_names)
    dollo = 1 + dollo_losses(mrca, presence_leaf_names)

    return {
        "is_monophyletic": bool(is_mono), "clade_type": clade_type,
        "clade_purity": clade_purity,
        "fitch_transitions": fitch, "dollo_transitions": dollo
    }


def write_itol_binary(path, clusters_to_plot, matrix, genome_to_leaf):
    """
    Write an iTOL DATASET_BINARY annotation file for the given clusters.
    Arguments:
        path (str): Path to the output file.
        clusters_to_plot (list): List of cluster IDs to include.
        matrix (dict): Presence/absence matrix.
        genome_to_leaf (dict): Mapping of genome IDs to tree leaf names.
    """
    field_labels = clusters_to_plot
    colors = [c for _, c in zip(field_labels, cycle(ITOL_PALETTE))]
    shapes = ["2"] * len(field_labels)  # 2 = filled circle

    with open(path, "w") as f:
        f.write("DATASET_BINARY\n")
        f.write("SEPARATOR TAB\n")
        f.write("DATASET_LABEL\tcluster_presence\n")
        f.write("COLOR\t#000000\n")
        f.write("FIELD_SHAPES\t" + "\t".join(shapes) + "\n")
        f.write("FIELD_LABELS\t" + "\t".join(field_labels) + "\n")
        f.write("FIELD_COLORS\t" + "\t".join(colors) + "\n")
        f.write("DATA\n")
        for genome, leaf_name in genome_to_leaf.items():
            if leaf_name is None:
                continue
            # For each cluster, get the presence/absence value for this genome (default to "-1" if not found)
            # replace 0 values with "-1" to indicate absence in iTOL
            values = [matrix[cid].get(genome, "-1") for cid in clusters_to_plot]
            values = ["-1" if v == "0" else v for v in values]
            f.write(leaf_name + "\t" + "\t".join(values) + "\n")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--presence-matrix", required=True,
                    help="presence_absence_matrix.tsv from analyse_cluster_spread.py")
    p.add_argument("--tree", required=True, help="Newick tree, leaf names = genome accessions")
    p.add_argument("--tree-metadata", required=True,
                    help="GTDB metadata file (gzip CSV) for mapping genome accessions to tree leaf names")
    p.add_argument("--output-report", required=True)
    p.add_argument("--output-itol", required=True)
    p.add_argument("--top-n", type=int, default=20,
                    help="Number of largest clusters to include in the iTOL annotation")
    p.add_argument("--expected-max-unmatched", type=int, default=100,
                    help="Expected number of matrix genomes with no matching tree leaf "
                         "(e.g. incomplete GenBank downloads). If the observed count "
                         "exceeds this, raise rather than continue silently, since it "
                         "likely indicates an accession-format mismatch (default: 100)")
    args = p.parse_args()

    genomes, matrix = load_presence_matrix(args.presence_matrix)
    tree = Tree(args.tree, format=1)
    log.info("Loaded matrix: %d clusters x %d genomes; tree: %d leaves",
              len(matrix), len(genomes), len(tree.get_leaves()))

    assembly_to_gtdb_dict = get_assembly_to_gtdb_dict(args.tree_metadata)

    leaves_to_keep, _ = ncbi_accessions_to_gtdb(genomes, assembly_to_gtdb_dict)
    tree = trim_tree_to_taxa(tree, leaves_to_keep)
    log.info("Pruned tree to %d leaves matching matrix genomes", len(tree.get_leaves()))

    tree_leaf_names = {leaf.name for leaf in tree.get_leaves()}
    genome_to_leaf, unmatched = build_genome_to_leaf(genomes, assembly_to_gtdb_dict, tree_leaf_names)

    if len(unmatched) > args.expected_max_unmatched:
        raise ValueError(
            f"{len(unmatched)} unmatched genomes exceeds --expected-max-unmatched "
            f"({args.expected_max_unmatched}). This is well above the expected dropout from "
            f"incomplete downloads - check for an accession-format mismatch between the "
            f"presence matrix and tree metadata before proceeding."
        )

    report_rows = []
    for cluster_id, genome_values in matrix.items():
        presence_leaf_names = {
            genome_to_leaf[g] for g, v in genome_values.items()
            if v == "1" and genome_to_leaf.get(g) is not None
        }
        stats = analyse_cluster(tree, presence_leaf_names)
        report_rows.append({
            "cluster_id": cluster_id,
            "n_genomes": len(presence_leaf_names),
            **stats,
        })

    report_rows.sort(key=lambda r: r["n_genomes"], reverse=True)

    fields = ["cluster_id", "n_genomes", "is_monophyletic", "clade_type",
              "clade_purity", "fitch_transitions", "dollo_transitions"]
    with open(args.output_report, "w") as f:
        f.write("\t".join(fields) + "\n")
        for r in report_rows:
            f.write("\t".join(str(r[fld]) for fld in fields) + "\n")

    top_clusters = [r["cluster_id"] for r in report_rows[:args.top_n]]
    write_itol_binary(args.output_itol, top_clusters, matrix, genome_to_leaf)

    log.info("Done. Wrote %s, %s", args.output_report, args.output_itol)


if __name__ == "__main__":
    main()

# example usage:
# python workflow/src/analyse_cluster_phylogeny.py \
#   --presence-matrix results/cluster_spread/presence_absence_matrix.tsv \
#   --tree results/phylogeny/gtdb_r207_bac120_r202 \
#   --tree-metadata data/gtdb_r207_metadata.csv.gz \
#   --output-report results/cluster_phylogeny/cluster_phylogeny_report.tsv \
#   --output-itol results/cluster_phylogeny/cluster_phylogeny_itol \
#   --top-n 20 \
#   --expected-max-unmatched 300
