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
from DarwinsRNAHunt.genome_analysis import parse_fasta, write_fasta_output
from DarwinsRNAHunt.sequence import ProteinSequence


def extract_cluster_proteins(clusters_file, cluster_id):
    """Extract protein IDs for a specific cluster from mmseqs2 TSV.
    Args:
        clusters_file: str filepath TSV format (mmseqs2 createtsv output):
            cluster_representative_id    cluster_member_id
        cluster_id: str
    Returns:
        set[str] ids of proteins that are members of the cluster
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

    output_path = Path(args.output)
    if proteins:
        protein_output = []
        for protein_id in proteins:
            record = sequences.get(protein_id)
            if record is None:
                print(f"Warning: {protein_id} not found in {args.fasta}, skipping")
                continue

            print(f"Extracting protein {protein_id}")
            # protein_id here is the full FASTA header (assembly|protein_id|type|
            # seqid:start-end(strand)), the same join-key format written by
            # ProteinSequence.to_header() in extract_flanking_cds.py -- so it can
            # be parsed directly rather than needing seqid/start/end/strand passed in.
            protein = ProteinSequence.from_header(protein_id, sequence=str(record.seq))
            if protein is None:
                print(f"Warning: could not parse header for {protein_id}, skipping")
                continue

            protein_output.append(protein)

        if protein_output:
            write_fasta_output(protein_output, output_path)
            print(f"Wrote {len(protein_output)} proteins to {args.output}")
        else:
            print(f"Warning: No proteins successfully extracted for cluster {args.cluster_id}")
            output_path.touch()
    else:
        print(f"Warning: No proteins found for cluster {args.cluster_id}")
        output_path.touch()


if __name__ == "__main__":
    main()