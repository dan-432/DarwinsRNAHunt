"""
author: Daniel Dachs
date: 16/02/2026
version: 1



"""

from interpro_access import *
from gtdb_access import *
import numpy as np
from copy import deepcopy
from ete3 import Tree
    

def trim_tree_to_taxa(tree, taxa_of_interest):
    """
    Prune GTDB tree to only include specified taxa.
    
    Args:
        tree: Bio.Phylo.Newick.Tree object (full GTDB tree)
        taxa_of_interest: list of GTDB IDs to keep (e.g., ['RS_GCF_950073225.1', ...])
        
    Returns:
        Bio.Phylo.Newick.Tree: Pruned tree containing only specified taxa
        
    """
    # Get all terminal nodes (species) in the tree
    all_terminals = tree.get_terminals()

    subtree = deepcopy(tree)
    
    # Remove unwanted terminals
    for terminal in all_terminals:
        if terminal.name not in taxa_of_interest:
            subtree.prune(terminal)
    
    # Clean up internal nodes that now have only one descendant
    # (Phylo doesn't do this automatically)
    #tree.collapse_all(lambda c: c.branch_length < 0.00001)
    
    return subtree

def calculate_distance_matrix(tree):
    """NOT WORKING
    Calculate pairwise distances between all taxa"""

    leaves = tree.get_terminals()

    leaf_names = [terminal.name for terminal in leaves]

    n = len(leaf_names)
    dist_matrix = np.zeros((n, n))
    
    for i, taxon1 in enumerate(leaves):
        for j, taxon2 in enumerate(leaves):
            if i != j:
                #print("{taxon1}, {taxon2}")
                dist_matrix[i, j] = tree.distance(taxon1, taxon2)
    
    return dist_matrix

def greedy_maxmin_selection(dist_matrix, taxa, n_select=100):
    """
    UNTESTED
    Greedy MaxMin algorithm to select maximally diverse subset
    
    At each step, add the taxon that maximizes the minimum distance
    to already-selected taxa
    """
    n_taxa = len(taxa)
    selected_indices = []
    
    # Start with a random taxon (or the most central one)
    selected_indices.append(0)
    
    # Greedily add taxa
    for _ in range(n_select - 1):
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
        
        selected_indices.append(best_candidate)
    
    return [taxa[i] for i in selected_indices]

def get_branch_length_to_root(tree, taxon_clade):
    """
    Get the branch length from root to a taxon.
    Fast: O(tree depth) per taxon
    """
    return tree.distance(tree.root, taxon_clade)


def get_branch_length_to_root(tree, taxon_clade):
    """
    Get the branch length from root to a taxon.
    Fast: O(tree depth) per taxon
    """
    return tree.distance(tree.root, taxon_clade)


def hierarchical_selection_from_tree(tree, taxa_names, n_select=100):
    """
    Select diverse taxa using tree structure directly - NO distance matrix!
    
    Strategy: Recursively split the tree and sample from each partition.
    Time: O(n log n) instead of O(n²)
    """
    print(f"Selecting {n_select} taxa using tree-based partitioning...")
    
    # Map names to clades
    name_to_clade = {t.name: t for t in tree.get_terminals()}
    taxa_clades = [name_to_clade[name] for name in taxa_names if name in name_to_clade]
    
    def partition_and_sample(clades, n_needed):
        """Recursively partition clades and sample evenly"""
        if n_needed >= len(clades):
            return clades
        
        if n_needed == 1:
            # Return most central clade
            root_dists = [tree.distance(tree.root, c) for c in clades]
            median_dist = np.median(root_dists)
            closest_to_median = min(range(len(clades)), 
                                   key=lambda i: abs(root_dists[i] - median_dist))
            return [clades[closest_to_median]]
        
        # Find the deepest common ancestor
        common_ancestor = tree.common_ancestor(clades)
        
        # Split into two groups based on which child clade they belong to
        if len(common_ancestor.clades) < 2:
            # Can't split further, random sample
            import random
            return random.sample(clades, n_needed)
        
        # Partition clades by which subtree they're in
        partitions = [[] for _ in common_ancestor.clades]
        for clade in clades:
            for i, child in enumerate(common_ancestor.clades):
                if child.is_parent_of(clade) or child == clade:
                    partitions[i].append(clade)
                    break
        
        # Remove empty partitions
        partitions = [p for p in partitions if p]
        
        # Sample proportionally from each partition
        selected = []
        samples_per_partition = [max(1, int(n_needed * len(p) / len(clades))) 
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
    
    selected_clades = partition_and_sample(taxa_clades, n_select)
    selected_names = [c.name for c in selected_clades]
    
    print(f"Selection complete! Selected {len(selected_names)} taxa")
    return selected_names



def main():
    # List of species taxids from the data we previouly downloaded
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


    #gtdb_tree = load_bacteria_tree(os.path.join(download_location, tree_file_name))

    # note this step removes non-bacterial taxa from our dataset as gtdb bacterial tree is just bacteria
    #bluf_tree = trim_tree_to_taxa(gtdb_tree, bluf_gtdb_taxaids)
    #save_tree(bluf_tree, "BLUF_bacGTDBtaxanomy.tree")

    bluf_tree = load_tree(os.path.join(download_location, "BLUF_bacGTDBtaxanomy.tree"))

    bluf_gtdb_taxaids = [terminal.name for terminal in bluf_tree.get_terminals()]
    print("loaded bluf tree from memory, calculating distance matrix between %d^2 taxonomies..."%len(bluf_gtdb_taxaids))

    #dist_matrix = calculate_distance_matrix(bluf_tree)
    #print("distance matrix calulated, doing selection...")
    #selected_100 = greedy_maxmin_selection(dist_matrix, bluf_gtdb_taxaids, n_select=100)

    #
    selected_100 = hierarchical_selection_from_tree(
        bluf_tree, 
        bluf_gtdb_taxaids, 
        n_select=100
    )

    print(selected_100)

    selected_100_tree = trim_tree_to_taxa(bluf_tree, selected_100)

    Phylo.draw(selected_100_tree)

    

if __name__ == "__main__":
    main()
