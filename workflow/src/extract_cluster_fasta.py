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
    
    if proteins:
        output_path = Path(args.output)
        protein_output = []
        for protein_id in proteins:
            seq = sequences[protein_id]

            # header eg: GCA_000026285.1|CAR04672.1|CDS|CU928161.2:3405150-3406199(+)
            assembly, prot_id, seq_type, seq_id, start, end, strand = re.split(r'|:()-', protein_id)
            protein_output.append({
                        "protein_id": f"{assembly}|{prot_id}",
                        "type": seq_type,
                        "seqid": seq_id,
                        "start": start,
                        "end": end,
                        "strand": strand,
                        "sequence": seq
                    })
        write_fasta_output(protein_output, output_path)
        print(f"Wrote {len(protein_output)} proteins to {args.output}")
    else:
        print(f"Warning: No proteins found for cluster {args.cluster_id}")
        Path(args.output).touch()
 
if __name__ == "__main__":
    main()