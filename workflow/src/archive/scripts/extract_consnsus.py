
"""
Extract consensus sequence from Stockholm alignment
"""
import argparse
from collections import Counter

def parse_stockholm(sto_file):
    """Parse Stockholm file and return sequences."""
    sequences = []
    with open(sto_file, 'r') as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('#') or line.startswith('//') or not line.strip():
                continue
            parts = line.split(None, 1)
            if len(parts) == 2:
                sequences.append(parts[1])
    return sequences

def get_consensus(sequences):
    """Get consensus sequence from alignment."""
    if not sequences:
        return ""
    
    consensus = []
    for i in range(len(sequences[0])):
        # Get column
        column = [seq[i] for seq in sequences if i < len(seq)]
        # Remove gaps
        bases = [b for b in column if b not in ['-', '.', '_']]
        if bases:
            # Most common base
            most_common = Counter(bases).most_common(1)[0][0]
            consensus.append(most_common)
        else:
            consensus.append('-')
    
    # Remove all gaps
    return ''.join(consensus).replace('-', '')

def main():
    parser = argparse.ArgumentParser(description='Extract consensus from Stockholm')
    parser.add_argument('--input-sto', required=True)
    parser.add_argument('--output-fasta', required=True)
    parser.add_argument('--motif-id', required=True)
    
    args = parser.parse_args()
    
    sequences = parse_stockholm(args.input_sto)
    consensus = get_consensus(sequences)
    
    with open(args.output_fasta, 'w') as f:
        f.write(f">{args.motif_id}_consensus\n")
        f.write(f"{consensus}\n")
    
    print(f"Consensus sequence ({len(consensus)} nt) written to {args.output_fasta}")

if __name__ == '__main__':
    main()