# workflow/src/remap_cluster_ids.py

import sys
import argparse

def main():
    parser = argparse.ArgumentParser(description="Remap complex cluster IDs to simple sequential IDs")
    parser.add_argument("--input", required=True, help="Input TSV (member_id, protein_member)")
    parser.add_argument("--output", required=True, help="Output TSV (sequential_cluster_id, protein_member)")
    args = parser.parse_args()
    
    seen = {}
    counter = 0

    print(f"Remapping cluster IDs from {args.input} to {args.output}", file=sys.stderr)
    
    with open(args.input) as infile, open(args.output, "w") as outfile:
        for line in infile:
            old_id, member = line.strip().split('\t')
            
            if old_id not in seen:
                seen[old_id] = counter
                counter += 1
            
            outfile.write(f"{seen[old_id]}\t{member}\n")

if __name__ == "__main__":
    main()