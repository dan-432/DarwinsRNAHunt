"""
Extract protein sequences for a specific cluster from mmseqs2 output.
 
Usage:
    python extract_cluster_fasta.py \
        --fasta combined_flanking_cds.faa \
        --clusters clusters.tsv \
        --cluster_id 0 \
        --output cluster_0.faa
"""
 
import argparse
from pathlib import Path
import re
from DarwinsRNAHunt.genome_analysis import parse_fasta, write_fasta_output
 
def extract_cluster_proteins(clusters_file, cluster_id):
    """Extract protein IDs for a specific cluster from mmseqs2 TSV.
    
    Args: 
        cluster_file str filepath TSV format (mmseqs2 createtsv output):
            cluster_representative_id    cluster_member_id
        cluster_id str
    Returns:
        [str] ids of proteins members of cluster
    """
    proteins = set()
    with open(clusters_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                representative = parts[0]
                member = parts[1]
                
                # Check if this line belongs to our cluster
                # mmseqs2 output has representative in first occurrence of cluster
                # We match by cluster ID (derived from representative)
                if representative == str(cluster_id):
                    proteins.add(member)
    
    return proteins

def parse_protein_id(protein_id, sequence):
    """Parse protein ID into its components.
    YUCKY GROSS REGEX (VOMIT) it works
    Example protein ID: GCA_000026285.1|CAR04672.1|CDS|CU928161.2:3405150-3406199(+)
    Returns:
        dict with keys: assembly, prot_id, seq_type, seq_id, start, end, strand
    """
    pattern = r'([^|]+)\|([^|]+)\|([^|]+)\|([^:]+):(\d+)-(\d+)\(([+-])\)'
    match = re.match(pattern, protein_id)
    if match:
        return {
            'protein_id': f"{match.group(1)}|{match.group(2)}",
            'type': match.group(3),
            'seqid': match.group(4),
            'start': int(match.group(5)),
            'end': int(match.group(6)),
            'strand': match.group(7),
            'sequence': sequence  # Fill in the actual sequence
        }
    return None

 
def main():
    parser = argparse.ArgumentParser(
        description="Extract proteins for a specific mmseqs2 cluster"
    )
    parser.add_argument("--fasta", required=True, help="Combined protein FASTA file")
    parser.add_argument("--clusters", required=True, help="mmseqs2 clusters TSV file")
    parser.add_argument("--cluster_id", required=True, help="Cluster ID to extract")
    parser.add_argument("--output", required=True, help="Output FASTA file")
    
    args = parser.parse_args()
    
    print(f"Parsing FASTA: {args.fasta}")
    sequences = parse_fasta(args.fasta)
    print(f"Loaded {len(sequences)} proteins")
    
    print(f"Extracting cluster {args.cluster_id} from {args.clusters}")
    proteins = extract_cluster_proteins(args.clusters, args.cluster_id)
    print(f"Found {len(proteins)} proteins in cluster {args.cluster_id}")

    output_path = Path(args.output)
    
    if proteins:
        protein_output = []
        for protein_id in proteins:
            print(f"Extracting protein {protein_id}")
            seq = str(sequences[protein_id].seq)

            protein_output.append(parse_protein_id(protein_id, seq))
        write_fasta_output(protein_output, output_path)
        print(f"Wrote {len(protein_output)} proteins to {args.output}")
    else:
        print(f"Warning: No proteins found for cluster {args.cluster_id}")
        output_path.touch()
 
if __name__ == "__main__":
    main()