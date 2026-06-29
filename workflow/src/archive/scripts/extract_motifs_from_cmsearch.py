import argparse
import pandas as pd

COLUMNS=["seq_id", "start", "end", "bit_score"]

def get_threshold(family_db, family_id):
    """Extract cmsearch threshold from family.txt for this family."""
    df = pd.read_table(family_db, header=None)
    row = df[df[0] == family_id]
    if row.empty:
        print("family not found using default bit score threshold")
        return 30.0  # fallback default
    return float(row.iloc[0, 1])

def parse_cmsearch_tblout(tblout_file, threshold):
    """Parse cmsearch --tblout output, filter by threshold."""
    motifs = []
    with open(tblout_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.split()
            if len(fields) < 16:
                continue
            seq_id = fields[0]
            bit_score = float(fields[14])
            start = int(fields[7])
            end = int(fields[8])
            
            if bit_score >= threshold:
                motifs.append((seq_id, start, end, bit_score))
    return motifs

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cmsearch_tblout", required=True)
    parser.add_argument("--family_db", required=True)
    parser.add_argument("--family", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    
    threshold = get_threshold(args.family_db, args.family)
    motifs = parse_cmsearch_tblout(args.cmsearch_tblout, threshold)
    
    df = pd.DataFrame(motifs, COLUMNS)
    df.to_csv(args.output, sep='\t', index=False)
    print(f"Found {len(motifs)} motifs for {args.family}")