from ete3 import Tree

# GTDB trees use quoted labels and sometimes non-standard bootstrap formatting,
# so format=1 (flexible) is usually more forgiving than the default.
tree = Tree("data/gtdb/bac120.tree", format=1, quoted_node_names=True)

# Get genome accessions for the family from the taxonomy TSV


family_tips = set()
with open("bac120_taxonomy-2.tsv") as f:
    for line in f:
        accession, lineage = line.rstrip("\n").split("\t")
        if "o__Lactobacillales" in lineage:
            family_tips.add(accession)

print(f"{len(family_tips)} genomes assigned to o__Lactobacillales")

# Match against actual tip names in the tree (accession formatting can differ
# slightly, e.g. GB_/RS_ prefixes) — check a few first if this comes up empty
tree_tips = set(tree.get_leaf_names())
matched = family_tips & tree_tips
print(f"{len(matched)} of those found as tips in the tree")

if matched:
    mrca = tree.get_common_ancestor(list(matched))
    print("MRCA node name:", mrca.name)
    print("Number of tips under MRCA:", len(mrca.get_leaves()))
    # Sanity check: are ALL genus/family tips under this MRCA, or did it
    # have to reach further up and pull in unrelated taxa (non-monophyletic)?
    extra = set(mrca.get_leaf_names()) - matched
    print("Non-family tips swept into this clade:", len(extra))
