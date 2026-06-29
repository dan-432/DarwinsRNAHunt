"""Proof of concept script, might become a module or speperate project later
For now the goal is to generate RNA motif sequences using a genetic algorithm and a covariance model. 

Logic:
1. Generate a starting sequence from the input covariance model (CM), using Infernal's cmemit or a custom implementation.
2. Create a population of sequences by mutating the starting sequence.
3. Evaluate the fitness of each sequence in the population using the covariance model, possibly by calculating the score of the sequence against the CM using Infernal's cmsearch or a custom implementation.
4. Select sequences based on their fitness scores to create a new generation of sequences.
5. Repeat the process for a specified number of generations or until a certain fitness threshold is reached.
6. Output the final set of generated RNA motif sequences in fasta format.
"""

import argparse
import math
import subprocess
import random
from Bio.Seq import Seq, MutableSeq
import json
from typing import Dict, List, Tuple, Optional

POSSIBLE_NUCLEOTIDES = ['A', 'U', 'C', 'G']

class RNA_Sequence:
    """Class to represent an RNA sequence."""
    def __init__(self, sequence):
        """Initialize the RNA sequence.
        Args:
            sequence (str or Seq or MutableSeq): The RNA sequence."""
        self.sequence = list(sequence)  # Store the sequence as a list of chars for easy manipulation
        self.parent_sequence = None  # To track the parent sequence for mutation history
        self.distance_from_parent = 0  # To track the number of mutations from the parent sequence
        self.fitness_score = None  # To store the fitness score of the sequence

    def mutate_to_child(self, mutation_rate):
        """Mutate the current sequence to create a child sequence.
        Args:
            mutation_rate (float): The rate at which mutations occur, between 0 and 1.
        Returns:
            RNA_Sequence: A new RNA_Sequence object representing the child sequence.
        """        
        child_sequence = self.sequence

        distance = 0

        for i in range(len(child_sequence)):
            if random.random() < mutation_rate:  # Mutate with the given mutation rate
                original_nucleotide = child_sequence[i]
                # Randomly choose a different nucleotide to mutate to
                possible_nucleotides = POSSIBLE_NUCLEOTIDES.copy()
                possible_nucleotides.remove(original_nucleotide)  # Remove the original nucleotide from the options
                child_sequence[i] = random.choice(possible_nucleotides)  # Mutate to a different nucleotide
                distance += 1

        child_rna_sequence = RNA_Sequence(child_sequence)  # Create a new RNA_Sequence object for the child
        child_rna_sequence.parent_sequence = self  # Set the parent sequence object (not the sequence itself)
        child_rna_sequence.distance_from_parent = distance
        return child_rna_sequence

    def set_fitness_score(self, score):
        """Set the fitness score of the sequence.
        Args:
            score (float): The fitness score to set."""
        self.fitness_score = score
    
    def get_fitness_score(self):
        """Get the fitness score of the sequence.
        Returns:
            float: The fitness score of the sequence."""
        return self.fitness_score

class RNA_Tree:
    """Class to represent a tree of RNA sequences for the genetic algorithm.
    This is an In-Tree, each sequence points to its parent up to the root sequence. Quicker more effiecient access to leaves and root."""
    def __init__(self, population_size_mu, population_size_sd, root_sequence):
        self.population_size_mu = population_size_mu
        self.population_size_sd = population_size_sd
        self.population = [RNA_Sequence(str(root_sequence))]  # Initialize the population with the root sequence

    def get_root_sequence(self):
        """Find and return the root sequence by traversing up from any leaf."""
        if not self.population:
            return None
        
        current = self.population[0]
        while current.parent_sequence is not None:
            current = current.parent_sequence
        
        return ''.join(current.sequence)

    def get_current_generation(self):
        """Return the current generation of sequences.
        Returns:
            list[Bio.Seq.Seq]: A list of RNA sequences in the population."""
        pop = []

        for sequence in self.population:
            # cast to Seq object
            pop.append(Seq(''.join(sequence.sequence)))

        return pop

    def mutate(self, mutation_rate):
        """Mutate the current population to create a new generation of sequences.
        Args:
            mutation_rate (float): The rate at which mutations occur, between 0 and 1.
        """
        next_generation = []

        for sequence in self.population:
            left_child = sequence.mutate_to_child(mutation_rate)
            right_child = sequence.mutate_to_child(mutation_rate)
            next_generation.append(left_child)
            next_generation.append(right_child)

        self.population = next_generation  # Update the population to the new generation

    def select(self, cm_file, scoring_function, selection_strength):
        """
        Evaluate the fitness of each sequence in the population using the covariance model.
        """

        # Calculate fitness scores for each sequence in the population using the provided scoring function and covariance model file
        for sequence in self.population:
            score = scoring_function(sequence, cm_file)
            sequence.set_fitness_score(score)

        # select sequences based on their fitness scores to create a new generation of sequences with probabalistic selection based on fitness scores and the selection strength parameter
        # this is an implementation of boltzmann selection, where the probability of selecting a sequence is proportional to its fitness score raised to the power of the selection strength parameter

        temperature = 1.0 / selection_strength if selection_strength > 0 else 1.0
        
        next_gen_size = int(random.gauss(self.population_size_mu, self.population_size_sd))

        if len(self.population) < next_gen_size:
            next_gen_size = len(self.population)  # Limit the next generation size to the current population size
            # JUST UNTILL WE GET TO A BG ENOUGH POPULATION. I DONT RLY LIKE THIS!!! GRRRRR
        next_generation = []
        
        fitness_scores = [seq.get_fitness_score() for seq in self.population]
        max_fitness = max(fitness_scores)
        
        # Boltzmann probabilities
        probabilities = [math.exp((score - max_fitness) / temperature) for score in fitness_scores]
        total_prob = sum(probabilities)
        probabilities = [p / total_prob for p in probabilities]
        
        for _ in range(next_gen_size):
            pick = random.random()
            current = 0
            for i, seq in enumerate(self.population):
                current += probabilities[i]
                if current >= pick:
                    next_generation.append(seq)
                    break

        self.population = next_generation  # Update the population to the new generation

    def visualize_tree_text(self, max_depth: Optional[int] = None) -> str:
        """
        Generate a text-based representation of the RNA tree.
        
        Args:
            max_depth: Maximum depth to display (None for all).
        
        Returns:
            String representation of the tree.
        """
        lines = []
        lines.append("RNA Tree Visualization (Text)")
        lines.append("=" * 80)
        
        # Create a mapping from parent sequence to children for quick lookup
        parent_to_children = {}
        seq_to_fitness = {}
        seq_id_map = {}
        counter = [0]  # Use list to allow modification in nested function
        
        # First pass: assign IDs and collect metadata
        def assign_ids(seq, visited=None):
            if visited is None:
                visited = set()
            seq_str = ''.join(seq.sequence)
            if seq_str in visited:
                return
            visited.add(seq_str)
            seq_id_map[seq_str] = counter[0]
            seq_to_fitness[seq_str] = seq.get_fitness_score()
            counter[0] += 1
        
        for seq in self.population:
            assign_ids(seq)
        
        # Build parent-to-children mapping by traversing tree
        for seq in self.population:
            seq_str = ''.join(seq.sequence)
            if seq.parent_sequence is not None:
                parent_str = ''.join(seq.parent_sequence.sequence)
                if parent_str not in parent_to_children:
                    parent_to_children[parent_str] = []
                parent_to_children[parent_str].append((seq_str, seq.distance_from_parent))
        
        # Tree traversal
        def print_tree(seq_str, prefix="", depth=0):
            if max_depth is not None and depth > max_depth:
                return
            
            seq_id = seq_id_map.get(seq_str, "?")
            fitness = seq_to_fitness.get(seq_str, None)
            fitness_str = f"fitness={fitness:.4f}" if fitness is not None else "fitness=?"
            
            lines.append(f"{prefix}├─ Seq#{seq_id}: {seq_str[:30]}... ({fitness_str})")
            
            if seq_str in parent_to_children:
                children = parent_to_children[seq_str]
                for i, (child_str, distance) in enumerate(children):
                    is_last = i == len(children) - 1
                    child_prefix = prefix + ("    " if is_last else "│   ")
                    lines.append(f"{prefix}│   └─[dist={distance}]")
                    print_tree(child_str, child_prefix, depth + 1)
        
        # Start from root
        root_str = ''.join(self.population[0].sequence)
        lines.append(f"Root: {root_str[:30]}... (fitness={self.population[0].get_fitness_score()})")
        
        # Print all sequences (since it's an inverted tree, we show all leaves)
        lines.append("\nCurrent Population:")
        for i, seq in enumerate(self.population):
            seq_str = ''.join(seq.sequence)
            distance = seq.distance_from_parent if seq.distance_from_parent else 0
            fitness = seq.get_fitness_score() if seq.get_fitness_score() is not None else 0
            parent_hint = f" [parent_dist={distance}]" if seq.parent_sequence else " [root]"
            lines.append(f"  {i}: {seq_str[:40]}...{parent_hint} (fitness={fitness:.4f})")
        
        return "\n".join(lines)

    def visualize_tree_newick(self, output_file: str = "rna_tree.nwk") -> str:
        """
        Generate a Newick format representation of the RNA tree by traversing parent chains.
        Each leaf in the population is traced back to the root, building a complete tree.
        
        Args:
            output_file: Path to save the Newick file.
        
        Returns:
            Newick format string with complete evolutionary tree.
        """
        # Build complete sequence map and parent-child relationships
        seq_to_id = {}
        seq_info = {}
        id_counter = [0]
        children_map = {}  # Maps parent_seq_str -> list of (child_seq_str, distance)
        parent_map = {}  # Maps child_seq_str -> parent_seq_str
        all_leaf_strs = set()  # Track which sequences are in the current population
        
        # Mark all current population members as leaves
        for leaf_seq in self.population:
            all_leaf_strs.add(''.join(leaf_seq.sequence))
        
        # Traverse from each population member up to root, collecting ALL sequences and relationships
        visited_pairs = set()  # Track which (parent, child) pairs we've seen
        
        for leaf_seq in self.population:
            current_seq = leaf_seq
            
            # Traverse ALL the way up to root without breaking
            while current_seq is not None:
                current_str = ''.join(current_seq.sequence)
                
                # Assign ID if needed
                if current_str not in seq_to_id:
                    seq_to_id[current_str] = id_counter[0]
                    id_counter[0] += 1
                    seq_info[current_str] = {
                        'fitness': current_seq.get_fitness_score(),
                        'distance': current_seq.distance_from_parent or 0,
                        'is_leaf': current_str in all_leaf_strs
                    }
                else:
                    # Mark as leaf if it's in the population
                    if current_str in all_leaf_strs:
                        seq_info[current_str]['is_leaf'] = True
                
                # Record parent-child relationship
                if current_seq.parent_sequence is not None:
                    parent_seq_obj = current_seq.parent_sequence
                    parent_str = ''.join(parent_seq_obj.sequence)
                    
                    pair = (parent_str, current_str)
                    if pair not in visited_pairs:
                        visited_pairs.add(pair)
                        parent_map[current_str] = parent_str
                        
                        if parent_str not in children_map:
                            children_map[parent_str] = []
                        children_map[parent_str].append((current_str, current_seq.distance_from_parent or 0))
                
                # Move to parent
                current_seq = current_seq.parent_sequence
        
        # Find root(s) - sequences with no parent
        roots = []
        for seq_str in seq_to_id.keys():
            if seq_str not in parent_map:
                roots.append(seq_str)
        
        # Build Newick recursively
        def build_newick(seq_str, visited=None):
            if visited is None:
                visited = set()
            
            if seq_str in visited:
                return ""
            visited.add(seq_str)
            
            seq_id = seq_to_id.get(seq_str, "?")
            distance = seq_info.get(seq_str, {}).get('distance', 0)
            
            # Check for children
            if seq_str in children_map:
                children_list = children_map[seq_str]
                child_trees = []
                for child_str, child_dist in children_list:
                    child_tree = build_newick(child_str, visited)
                    if child_tree:
                        child_trees.append(child_tree)
                
                if child_trees:
                    children_str = ",".join(child_trees)
                    return f"({children_str})seq{seq_id}:{distance}"
            
            # Leaf node
            return f"seq{seq_id}:{distance}"
        
        # Build complete tree
        if len(roots) == 0:
            newick_str = "()empty:0;"
        elif len(roots) == 1:
            tree_str = build_newick(roots[0])
            newick_str = tree_str + ";"
        else:
            # Multiple roots - create super-root
            root_trees = []
            for root in roots:
                tree = build_newick(root)
                if tree:
                    root_trees.append(tree)
            newick_str = "(" + ",".join(root_trees) + ")superroot:0;"
        
        # Create header
        header = "# Newick format RNA evolutionary tree\n"
        header += "# Built by traversing parent chains from current population to root\n"
        header += "# Branch lengths represent mutation distances (Hamming distance from parent)\n"
        header += "# Leaf nodes (marked with *) are current population members\n"
        header += f"# Total nodes in tree: {len(seq_to_id)}\n"
        header += f"# Current population size (leaves): {len([s for s in seq_info.values() if s['is_leaf']])}\n"
        header += f"# Number of roots: {len(roots)}\n"
        header += "# \n# Node information (seq_id: fitness, parent_dist, is_leaf):\n"
        
        for seq_str, seq_id in sorted(seq_to_id.items(), key=lambda x: x[1]):
            fitness = seq_info[seq_str]['fitness']
            distance = seq_info[seq_str]['distance']
            is_leaf = seq_info[seq_str]['is_leaf']
            fitness_str = f"{fitness:.4f}" if fitness is not None else "?"
            leaf_marker = "*" if is_leaf else " "
            header += f"# seq{seq_id}: fitness={fitness_str}, dist={distance:2d} {leaf_marker}\n"
        
        header += "\n"
        output = header + newick_str
        
        # Write to file
        with open(output_file, 'w') as f:
            f.write(output)
        
        return output

    def export_tree_json(self, output_file: str = "rna_tree.json") -> Dict:
        """
        Export tree structure and metrics to JSON format.
        
        Args:
            output_file: Path to save the JSON file.
        
        Returns:
            Dictionary containing tree data.
        """
        tree_data = {
            "population_size_mu": self.population_size_mu,
            "population_size_sd": self.population_size_sd,
            "root_sequence": self.get_root_sequence(),
            "current_population_size": len(self.population),
            "sequences": []
        }
        
        for i, seq in enumerate(self.population):
            seq_info = {
                "id": i,
                "sequence": ''.join(seq.sequence),
                "fitness_score": seq.get_fitness_score(),
                "distance_from_parent": seq.distance_from_parent,
                "parent_sequence": ''.join(seq.parent_sequence.sequence) if seq.parent_sequence else None
            }
            tree_data["sequences"].append(seq_info)
        
        # Calculate statistics
        fitness_scores = [seq.get_fitness_score() for seq in self.population if seq.get_fitness_score() is not None]
        if fitness_scores:
            tree_data["fitness_stats"] = {
                "min": min(fitness_scores),
                "max": max(fitness_scores),
                "mean": sum(fitness_scores) / len(fitness_scores),
                "count": len(fitness_scores)
            }
        
        # Write to file
        with open(output_file, 'w') as f:
            json.dump(tree_data, f, indent=2)
        
        return tree_data

    def get_tree_statistics(self) -> Dict:
        """
        Calculate statistics about the current tree state.
        
        Returns:
            Dictionary containing tree statistics.
        """
        fitness_scores = [seq.get_fitness_score() for seq in self.population if seq.get_fitness_score() is not None]
        distances = [seq.distance_from_parent for seq in self.population if seq.distance_from_parent is not None]
        
        stats = {
            "population_size": len(self.population),
            "fitness": {
                "min": min(fitness_scores) if fitness_scores else None,
                "max": max(fitness_scores) if fitness_scores else None,
                "mean": sum(fitness_scores) / len(fitness_scores) if fitness_scores else None,
                "count": len(fitness_scores)
            },
            "mutations": {
                "min_distance": min(distances) if distances else None,
                "max_distance": max(distances) if distances else None,
                "mean_distance": sum(distances) / len(distances) if distances else None
            }
        }
        return stats

def generate_starting_sequence(cm_file, score_threshold):
    # Implement logic to generate a starting sequence from the covariance model using Infernal's cmemit or a custom implementation
    pass

def score_sequence_infernal_cmsearch(sequence, cm_file):
    # Implement logic to score a sequence against the covariance model using Infernal's cmsearch or a custom implementation
    pass

def score_sequence_test(sequence, cm_file):
    # Placeholder scoring function for testing purposes, returns a random score
    # high score closest to 50% GC content, low score furthest from 50% GC content

    return 1 - abs(gc_content(sequence) - 50) / 50  # Normalize to a score between 0 and 1, where 1 is best at 50% GC content

def gc_content(sequence):
    """Calculate the GC content of a given RNA sequence.
    Args:
        sequence (str or Seq): The RNA sequence for which to calculate GC content.
    Returns:
        float: The GC content of the sequence as a percentage."""
    gc_count = sequence.sequence.count('G') + sequence.sequence.count('C')
    total_count = len(sequence.sequence)
    return (gc_count / total_count) * 100 if total_count > 0 else 0


def main():

    parser = argparse.ArgumentParser(description="Generate RNA motif sequences using genetic algorithm")
    parser.add_argument("--cm_file", type=str, help="Path to covariant model file to build from")
    parser.add_argument("--population_size", type=int, default=100, help="Size of the population for the genetic algorithm")
    parser.add_argument("--generations", type=int, default=100, help="Number of generations to run the genetic algorithm")
    parser.add_argument("--mutation_rate", type=float, default=0.01, help="Mutation rate for the genetic algorithm")
    parser.add_argument("--selection_strength", type=float, default=0.7, help="Strength of selection for the genetic algorithm, between 0 and 1")
    parser.add_argument("--output_file", type=str, help="File to save the generated RNA, fasta format")

    args = parser.parse_args()

    

if __name__ == "__main__":
    # test 

    rtree = RNA_Tree(100, 10, "AUGCUAGCUAGCCUGACUGAGCGCUAUAUGCGCAUAUGCGCAUGAGCUAGUCGUGCUAGAUGUCGCAUGAUGCUCGUAGCGAGCUAUGCCGACUCUGAGACGACUAGCGGCUAUAGCGC")

    for i in range(50):
        rtree.mutate(0.1)
        rtree.select("test.cm", score_sequence_test, 0.7)
        
        
    # Print text visualization
    print(f"\n{rtree.visualize_tree_text()}")

    rtree.visualize_tree_newick(f"rna_tree_gen_{i+1}.nwk")
    rtree.export_tree_json(f"rna_tree_gen_{i+1}.json")
    print(f"\nExported tree visualization to rna_tree_gen_{i+1}.nwk and rna_tree_gen_{i+1}.json")