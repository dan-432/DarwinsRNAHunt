
"""
Summarize Rfam annotation results
"""
import argparse
from collections import Counter

def parse_rfam_hits(tblout_file):
    """Parse cmscan tblout file."""
    hits = []
    with open(tblout_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 4:
                hits.append({
                    'target': parts[0],
                    'query': parts[1],
                    'accession': parts[2],
                    'evalue': parts[3],
                    'score': parts[4] if len(parts) > 4 else 'N/A'
                })
    return hits

def main():
    parser = argparse.ArgumentParser(description='Summarize Rfam hits')
    parser.add_argument('--rfam-hits', required=True)
    parser.add_argument('--motif-id', required=True)
    parser.add_argument('--seq-count', required=True)
    parser.add_argument('--output-summary', required=True)
    parser.add_argument('--output-flag', required=True)
    
    args = parser.parse_args()
    
    hits = parse_rfam_hits(args.rfam_hits)
    
    # Count unique families
    families = [h['query'] for h in hits]
    family_counts = Counter(families)
    unique_families = len(family_counts)
    total_hits = len(hits)
    
    # Write summary
    with open(args.output_summary, 'w') as f:
        f.write("Rfam Annotation Summary\n")
        f.write("======================\n")
        f.write(f"Motif: {args.motif_id}\n")
        f.write(f"Sequences searched: {args.seq_count}\n")
        f.write(f"\n")
        f.write(f"Total hits: {total_hits}\n")
        f.write(f"Unique Rfam families: {unique_families}\n")
        f.write(f"\n")
        
        if unique_families > 0:
            f.write("Status: KNOWN MOTIF\n")
            f.write("\n")
            f.write("Matched Rfam families:\n")
            for family, count in family_counts.most_common(10):
                f.write(f"  {count:4d}  {family}\n")
            f.write("\n")
            f.write("Top hits:\n")
            for hit in hits[:10]:
                f.write(f"  {hit['target']:15s} {hit['query']:20s} "
                       f"E={hit['evalue']:12s} Score={hit['score']:8s}\n")
        else:
            f.write("Status: NOVEL MOTIF\n")
            f.write("\n")
            f.write("No significant matches to known Rfam families.\n")
            f.write("This motif may represent a novel RNA structure.\n")
    
    # Write flag file
    with open(args.output_flag, 'w') as f:
        if unique_families > 0:
            f.write(f"KNOWN - matches {unique_families} Rfam families\n")
        else:
            f.write("NOVEL - no Rfam matches\n")
    
    print(f"Summary written to {args.output_summary}")
    print(f"Status: {'KNOWN' if unique_families > 0 else 'NOVEL'} "
          f"({unique_families} Rfam families, {total_hits} hits)")

if __name__ == '__main__':
    main()