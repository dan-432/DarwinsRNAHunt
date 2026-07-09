"""
author: Daniel Dachs
date: 16/02/2026
version: 2

Now using ete3 rather than bio phyla
"""

import numpy as np
from copy import deepcopy
from ete3 import Tree


def trim_tree_to_taxa(tree, taxa_of_interest):
    """
    Prune tree to only include specified taxa.
    Args:
        tree: ete3.Tree object
        taxa_of_interest: list of GTDB IDs to keep (e.g., ['RS_GCF_950073225.1', ...])
        representative_map: optional dict mapping GTDB ID -> GTDB species representative ID.
            If a taxon is missing from the tree and this is supplied, the taxon's
            representative will be substituted in instead (if the representative
            is itself in the tree).
    Returns:
        ete3.Tree: Pruned tree containing only specified taxa (or their representatives)
    """
    print(f"Trimming tree to {len(taxa_of_interest)} taxa...")

    subtree = deepcopy(tree)
    tree_leaves = {leaf.name for leaf in tree.iter_leaves()}

    taxa_in_tree = set()
    substituted = []
    taxa_not_in_tree = []

    for taxon in taxa_of_interest:
        if taxon in tree_leaves:
            taxa_in_tree.add(taxon)
            continue

        taxa_not_in_tree.append(taxon)

    print(f"  Taxa in tree: {len(taxa_in_tree) - len(substituted)}/{len(taxa_of_interest)} directly")

    if substituted:
        print(f"  Substituted {len(substituted)} missing taxa with their GTDB representative:")
        for taxon, rep in substituted[:10]:
            print(f"    - {taxon} -> {rep}")
        if len(substituted) > 10:
            print(f"    ... and {len(substituted) - 10} more")

    if taxa_not_in_tree:
        print(f"  WARNING: {len(taxa_not_in_tree)} taxa not found in tree (no usable representative)")
        for taxon in taxa_not_in_tree[:10]:
            print(f"    - {taxon}")
        if len(taxa_not_in_tree) > 10:
            print(f"    ... and {len(taxa_not_in_tree) - 10} more")

    if not taxa_in_tree:
        raise ValueError("ERROR: None of the specified taxa (or their representatives) were found in the tree!")

    subtree.prune(list(taxa_in_tree), preserve_branch_length=True)
    print(f"Trimmed tree has {len(list(subtree.iter_leaves()))} leaves")
    return subtree

def trim_tree_to_order(tree, taxorder_of_interest):
    """
    Prune tree to only include specified taxanomic order.
    
    Args:
        tree: ete3.Tree object
        taxorder_of_interest: taxanomic order of interest (GTDB nomenclature) to trim down to
        
    Returns:
        ete3.Tree: Subtree of specified taxanomic order
    """

    subtree = None

    for node in tree.traverse():
        if node.name == taxorder_of_interest:
            subtree = node
            break

    if not subtree:
        raise Exception(f"{taxorder_of_interest} could not be found in tree :(")

    return subtree

def get_tree_accession_ids(tree):
    """
    Returns genome accession ids of tree.
    """
    return tree.get_leaf_names()


def calculate_distance_matrix(tree):
    """
    Calculate pairwise distances between all taxa.
    
    Args:
        tree: ete3.Tree object
        
    Returns:
        np.ndarray: Distance matrix
    """
    leaves = list(tree.iter_leaves())
    leaf_names = [leaf.name for leaf in leaves]
    
    n = len(leaves)
    dist_matrix = np.zeros((n, n))
    
    print(f"Calculating {n}x{n} distance matrix...")
    
    for i, leaf1 in enumerate(leaves):
        if i % 100 == 0:
            print(f"  Progress: {i}/{n}")
        for j, leaf2 in enumerate(leaves):
            if i != j:
                # ete3's get_distance calculates phylogenetic distance
                dist_matrix[i, j] = tree.get_distance(leaf1, leaf2)
    
    print("Distance matrix complete!")
    return dist_matrix


def greedy_maxmin_selection(dist_matrix, taxa_names, n_select=100):
    """
    Greedy MaxMin algorithm to select maximally diverse subset.
    
    At each step, add the taxon that maximizes the minimum distance
    to already-selected taxa.
    
    Args:
        dist_matrix: np.ndarray of pairwise distances
        taxa_names: list of taxon names
        n_select: number to select
        
    Returns:
        list: Selected taxon names
    """
    n_taxa = len(taxa_names)
    selected_indices = []
    
    # Start with most central taxon
    central_idx = dist_matrix.sum(axis=1).argmin()
    selected_indices.append(central_idx)
    
    print(f"Selecting {n_select} taxa using MaxMin algorithm...")
    
    # Greedily add taxa
    for iteration in range(n_select - 1):
        if iteration % 10 == 0:
            print(f"  Progress: {iteration + 1}/{n_select}")
        
        max_min_dist = -1
        best_candidate = None
        
        for candidate in range(n_taxa):
            if candidate in selected_indices:
                continue
            
            # Find minimum distance from candidate to all selected taxa
            min_dist = min(dist_matrix[candidate, selected] 
                          for selected in selected_indices)
            
            # Keep track of candidate with maximum min-distance
            if min_dist > max_min_dist:
                max_min_dist = min_dist
                best_candidate = candidate
        
        if best_candidate is not None:
            selected_indices.append(best_candidate)
    
    return [taxa_names[i] for i in selected_indices]


def get_branch_length_to_root(tree, leaf):
    """
    Get the branch length from root to a leaf.
    
    Args:
        tree: ete3.Tree object
        leaf: ete3.TreeNode (leaf node)
        
    Returns:
        float: Distance from root to leaf
    """
    return tree.get_distance(tree, leaf)


def hierarchical_selection_from_tree(tree, taxa_names, n_select=100):
    """
    Select diverse taxa using tree structure directly - NO distance matrix!
    
    Strategy: Recursively split the tree and sample from each partition.
    Time: O(n log n) instead of O(n²)
    
    Args:
        tree: ete3.Tree object
        taxa_names: list of taxon names to select from
        n_select: number to select
        
    Returns:
        list: Selected taxon names
    """
    print(f"Selecting {n_select} taxa using tree-based partitioning...")
    
    # Map names to leaf nodes
    name_to_leaf = {leaf.name: leaf for leaf in tree.iter_leaves()}
    taxa_leaves = [name_to_leaf[name] for name in taxa_names if name in name_to_leaf]
    
    if len(taxa_leaves) == 0:
        print("ERROR: No taxa found in tree")
        return []
    
    def partition_and_sample(leaves, n_needed):
        """Recursively partition leaves and sample evenly"""
        if n_needed >= len(leaves):
            return leaves
        
        if n_needed == 1:
            # Return most central leaf (closest to median distance from root)
            root_dists = [tree.get_distance(tree, leaf) for leaf in leaves]
            median_dist = np.median(root_dists)
            closest_to_median = min(range(len(leaves)), 
                                   key=lambda i: abs(root_dists[i] - median_dist))
            return [leaves[closest_to_median]]
        
        # Find the common ancestor of all these leaves
        common_ancestor = tree.get_common_ancestor(leaves)
        
        # Get children of common ancestor
        children = common_ancestor.get_children()
        
        if len(children) < 2:
            # Can't split further, random sample
            import random
            random.seed(42)
            return random.sample(leaves, n_needed)
        
        # Partition leaves by which child subtree they belong to
        partitions = [[] for _ in children]
        for leaf in leaves:
            for i, child in enumerate(children):
                # Check if this leaf descends from this child
                if leaf in child.iter_leaves():
                    partitions[i].append(leaf)
                    break
        
        # Remove empty partitions
        partitions = [p for p in partitions if p]
        
        if len(partitions) == 0:
            import random
            random.seed(42)
            return random.sample(leaves, n_needed)
        
        # Sample proportionally from each partition
        selected = []
        samples_per_partition = [max(1, int(n_needed * len(p) / len(leaves))) 
                                for p in partitions]
        
        # Adjust to ensure we get exactly n_needed
        while sum(samples_per_partition) < n_needed:
            # Add to largest partition
            largest = max(range(len(partitions)), key=lambda i: len(partitions[i]))
            samples_per_partition[largest] += 1
        
        while sum(samples_per_partition) > n_needed:
            # Remove from partition with most samples
            largest = max(range(len(samples_per_partition)), 
                         key=lambda i: samples_per_partition[i])
            samples_per_partition[largest] -= 1
        
        # Recursively sample from each partition
        for partition, n_from_partition in zip(partitions, samples_per_partition):
            if n_from_partition > 0:
                selected.extend(partition_and_sample(partition, n_from_partition))
        
        return selected
    
    selected_leaves = partition_and_sample(taxa_leaves, n_select)
    selected_names = [leaf.name for leaf in selected_leaves]
    
    print(f"Selection complete! Selected {len(selected_names)} taxa")
    return selected_names


def test():

    from interpro_access import load_taxonomic_info, get_taxa_ids
    from gtdb_access import get_ncbi_to_gtdb_dict, load_tree, save_tree

    # List of species taxids from the data we previously downloaded
    bluf_taxa_info = load_taxonomic_info()
    bluf_ncbi_taxids = get_taxa_ids(bluf_taxa_info)
    
    taxa_dict = get_ncbi_to_gtdb_dict()
    
    bluf_gtdb_taxaids = []
    no_mapping_count = 0
    
    for ncbi_id in bluf_ncbi_taxids:
        if ncbi_id in taxa_dict.keys():
            bluf_gtdb_taxaids.extend(taxa_dict[ncbi_id])
        else:
            no_mapping_count += 1
    
    print(f"NCBI taxonomy IDs which do not map to GTDB accessions: {no_mapping_count} of a total {len(bluf_ncbi_taxids)}")
    
    # Option 1: Load full tree and trim (first time)
    # gtdb_tree = load_tree_ete3(os.path.join(download_location, tree_file_name))
    # bluf_tree = trim_tree_to_taxa(gtdb_tree, bluf_gtdb_taxaids)
    # save_tree_ete3(bluf_tree, os.path.join(download_location, "BLUF_bacGTDBtaxanomy.tree"))
    
    # Option 2: Load pre-saved trimmed tree
    bluf_tree = load_tree(os.path.join("data/trees", "BLUF_bacGTDBtaxanomy.tree"))
    
    # Get all taxa names from tree
    bluf_gtdb_taxaids = [leaf.name for leaf in bluf_tree.iter_leaves()]
    print(f"Loaded BLUF tree with {len(bluf_gtdb_taxaids)} taxa")
    
    # Select 100 diverse taxa
    selected_100 = hierarchical_selection_from_tree(
        bluf_tree, 
        bluf_gtdb_taxaids, 
        n_select=100
    )
    
    print(f"\nSelected {len(selected_100)} taxa:")
    for i, taxon in enumerate(selected_100[:10]):  # Show first 10
        print(f"  {i+1}. {taxon}")
    print(f"  ... and {len(selected_100) - 10} more")
    
    # Create tree with just selected taxa
    selected_100_tree = trim_tree_to_taxa(bluf_tree, selected_100)
    
    # Save selected tree
    save_tree(selected_100_tree, "data/trees/BLUF_selected_100.tree")
    
    # Visualize using ete3
    print("\nRendering tree...")
    selected_100_tree.show()  # Interactive viewer
    
    # Or save as image
    from ete3 import TreeStyle
    ts = TreeStyle()
    ts.show_branch_length = True
    ts.show_leaf_name = True
    selected_100_tree.render("data/selected_100_tree.pdf", tree_style=ts)
    print("Tree saved to data/selected_100_tree.pdf")


if __name__ == "__main__":
    test()