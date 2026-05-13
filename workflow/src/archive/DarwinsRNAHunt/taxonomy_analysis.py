"""
author: Daniel Dachs
date: 13/02/2026
version: 1

Huerta-Cepas, J., et al. (2016). "ETE 3: Reconstruction, Analysis, and Visualization of Phylogenomic Data." Molecular biology and evolution 33.

"""

from ete3 import NCBITaxa, NodeStyle, TextFace
from interpro_access import *
from collections import Counter
import matplotlib.pyplot as plt

def name_tree(tree, ncbi):
    for node in tree.traverse():
    
        node.name = ncbi.get_taxid_translator([int(node.name)])[int(node.name)]

    return tree

def get_family_tree(tree, ncbi):
    families = []

    for node in tree.traverse():
        rank = ncbi.get_rank([node.taxid])[node.taxid]
        if rank == "family":
            families.append(node.taxid)

    family_tree = ncbi.get_topology(families)

    family_tree = name_tree(family_tree, ncbi)

    return family_tree

def greedy_pd_ranking(tree, n_select=300):
    """
    Implementation of greedy phylogenetic diversity maximization.
    
    Iteratively selects the leaf whose addition contributes the most
    new branch length (uncovered by already-selected leaves).
    
    Returns a ranked list of (leaf_name, marginal_pd_contribution) tuples.
    """

    # Track which nodes have been "covered" (their branch length already counted)
    covered_nodes = set()
    selected = []
    remaining_leaves = set(tree.get_leaves())

    def marginal_pd(leaf):
        """Branch length contributed by adding this leaf to the current selection."""
        node = leaf
        total = 0.0
        while node is not None:
            if node in covered_nodes:
                break  # Everything above here is already counted
            total += node.dist
            node = node.up
        return total

    # Seed with the leaf on the longest root-to-tip path (maximizes first pick)
    first_leaf = max(remaining_leaves, key=lambda l: sum(n.dist for n in l.iter_ancestors()) + l.dist)
    
    # Mark its path to root as covered
    def cover_path(leaf):
        node = leaf
        while node is not None:
            if node in covered_nodes:
                break
            covered_nodes.add(node)
            node = node.up

    cover_path(first_leaf)
    selected.append((first_leaf.name, marginal_pd(first_leaf)))  # contribution = full path length
    # Recalculate first leaf's actual contribution
    selected[0] = (first_leaf.name, sum(n.dist for n in first_leaf.iter_ancestors()) + first_leaf.dist)
    remaining_leaves.remove(first_leaf)

    # Greedy iterations
    while len(selected) < n_select and remaining_leaves:
        best_leaf = max(remaining_leaves, key=marginal_pd)
        contribution = marginal_pd(best_leaf)
        
        if contribution == 0:
            # All remaining leaves are on fully-covered branches (shouldn't happen with a real tree)
            break
        
        cover_path(best_leaf)
        selected.append((best_leaf.name, contribution))
        remaining_leaves.remove(best_leaf)

        if len(selected) % 50 == 0:
            print(f"Selected {len(selected)} taxa, last contribution: {contribution:.4f}")

    return selected

def graph_PD_distribution(scores, show = False):

    # See where the natural breaks are
    score_counts = Counter(round(s, 1) for s in scores)
    print(sorted(score_counts.items(), reverse=True))

    # Plot the "scree plot" - look for the elbow
    plt.plot(range(len(scores)), scores)
    plt.xlabel("Rank")
    plt.ylabel("Marginal PD contribution")
    plt.axhline(y=2, color='r', linestyle='--', label='cutoff at 2')
    
    plt.savefig("output/PDRankingBLUF.png")
    if show: plt.show()

    # Cumulative PD captured - useful for justifying your cutoff
    cumulative = [sum(scores[:i]) for i in range(1, len(scores)+1)]
    total_pd = sum(s for _, s in greedy_pd_ranking(tree, n_select=len(tree.get_leaves())))
    plt.plot([c/total_pd * 100 for c in cumulative])
    plt.xlabel("Number of taxa selected")
    plt.ylabel("% total PD captured")
    
    plt.savefig("output/PDCoverageBLUF.png")
    if show: plt.show()


if __name__ == "__main__":

    #make NCBITaxa instance
    ncbi = NCBITaxa()

    # List of species taxids from the data we previouly downloaded
    taxa_info = load_taxonomic_info()

    taxids = get_taxa_ids(taxa_info)

    # Get the tree topology
    tree = ncbi.get_topology(taxids)

    tree = name_tree(tree, ncbi)

    # include only bactria
    bacteria_tree = None
    for node in tree.traverse():
        # Check if '2' (Bacteria) is in the lineage of this node
        if 2 in node.lineage:
            bacteria_tree = node
            break


    ranked = greedy_pd_ranking(bacteria_tree, n_select=200)

    high_ranked_node_names, scores = zip(*ranked)


    print(high_ranked_node_names)

    graph_PD_distribution(scores)

    # prune tree for high ranked nodes
    filtered_tree = tree.copy()

    filtered_tree.prune(high_ranked_node_names)

    # Style for internal/other nodes
    nstyle = NodeStyle()
    nstyle["fgcolor"] = "red"
    nstyle["size"] = 10

    for n in tree.traverse():
        if n.name in high_ranked_node_names:
            n.set_style(nstyle)
            # Optional: Add text face for names, as default names might be disabled
            name_face = TextFace(n.name, fgcolor="blue", fsize=12)
            n.add_face(name_face, column=0, position='branch-right')
    
    tree.render("output/SelectedSpeciesByPD.png")


    print(filtered_tree)