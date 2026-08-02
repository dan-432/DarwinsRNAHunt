"""
Extract nucleotide sequences for a specific cluster from mmseqs2 output.
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
from DarwinsRNAHunt.sequence import NucleotideSequence

def extract_cluster_sequences(clusters_file, cluster_id):
    """Extract sequence IDs for a specific cluster from mmseqs2 TSV.
    Args:
        clusters_file: str filepath TSV format (mmseqs2 createtsv output):
            cluster_representative_id    cluster_member_id
        cluster_id: str
    Returns:
        set[str] ids of sequences that are members of the cluster
    """
    sequences = set()
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
                    sequences.add(member)
    return sequences


def main():
    parser = argparse.ArgumentParser(
        description="Extract sequences for a specific mmseqs2 cluster"
    )
    parser.add_argument("--fasta", required=True, help="Combined sequence FASTA file")
    parser.add_argument("--clusters", required=True, help="mmseqs2 clusters TSV file")
    parser.add_argument("--cluster_id", required=True, help="Cluster ID to extract")
    parser.add_argument("--output", required=True, help="Output FASTA file")
    args = parser.parse_args()

    print(f"Parsing FASTA: {args.fasta}")
    full_sequences = parse_fasta(args.fasta)
    print(f"Loaded {len(full_sequences)} sequences")

    print(f"Extracting cluster {args.cluster_id} from {args.clusters}")
    sequences = extract_cluster_sequences(args.clusters, args.cluster_id)
    print(f"Found {len(sequences)} sequences in cluster {args.cluster_id}")

    output_path = Path(args.output)
    if sequences:
        sequence_output = []
        for sequence_id in sequences:
            record = full_sequences.get(sequence_id)
            if record is None:
                print(f"Warning: {sequence_id} not found in {args.fasta}, skipping")
                continue

            print(f"Extracting sequence {sequence_id}")
            # sequence_id here is the full FASTA header (assembly|sequence_id|type|
            # seqid:start-end(strand)), the same join-key format written by
            # NucleotideSequence.to_header() in extract_flanking_cds.py -- so it can
            # be parsed directly rather than needing seqid/start/end/strand passed in.
            sequence = NucleotideSequence.from_header(sequence_id, sequence=str(record.seq))
            if sequence is None:
                print(f"Warning: could not parse header for {sequence_id}, skipping")
                continue

            sequence_output.append(sequence)

        if sequence_output:
            write_fasta_output(sequence_output, output_path)
            print(f"Wrote {len(sequence_output)} sequences to {args.output}")
        else:
            print(f"Warning: No sequences successfully extracted for cluster {args.cluster_id}")
            output_path.touch()
    else:
        print(f"Warning: No sequences found for cluster {args.cluster_id}")
        output_path.touch()


if __name__ == "__main__":
    main()