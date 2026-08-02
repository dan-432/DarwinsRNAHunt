"""
Extract motif sequences from a FASTA file, filtered to the sequence IDs
that had motif hits in a mmseqs tblout file.
"""

import argparse
import json
from DarwinsRNAHunt.genome_analysis import parse_tblout, parse_fasta, write_fasta_output
from DarwinsRNAHunt.sequence import NucleotideSequence


def main():
    # example usage:
    # python workflow/src/extract_motif_sequences.py --fasta-file results/02_sequence_selection/03_combined/{target_domain}_domain_flanking_upstream.fasta --motif-hits-tbl results/03_motif_discovery/{target_domain}/02_motif_hits/{motif_id}_hits.tblout --output-fasta results/03_motif_discovery/{target_domain}/03_homolog_search/expanded_seqs/{motif_id}_expanded.fasta --output-ids results/03_motif_discovery/{target_domain}/03_homolog_search/expanded_seqs/{motif_id}_ids.json
    parser = argparse.ArgumentParser(description="Extract motif sequences from a FASTA file.")
    parser.add_argument("--fasta-file", help="Path to the input FASTA file.")
    parser.add_argument("--motif-hits-tbl", help="Path to the motif hits table.")
    parser.add_argument("--output-fasta", help="Path to the output FASTA file.")
    parser.add_argument("--output-ids", help="Path to the output file for sequence IDs.")
    args = parser.parse_args()

    print(f"Extracting motif sequences from {args.fasta_file} based on hits in {args.motif_hits_tbl}...")
    seqs = parse_fasta(args.fasta_file)
    motif_hits = parse_tblout(args.motif_hits_tbl)

    print(f"Found {len(motif_hits)} motif hits. Extracting corresponding sequences...")

    filtered_records = []

    for seq_id in motif_hits:
        if seq_id in seqs:
            record = seqs[seq_id]
            # record.id is expected to be a header in
            # 'protein_id|type|seqid:start-end(strand)' format, as written by
            # NucleotideSequence.to_header() upstream -- from_header parses
            # it and dispatches to the right subclass (CDSSequence/UTR5Sequence/
            # UTR3Sequence/Motif) based on the type field.
            parsed = NucleotideSequence.from_header(record.id, sequence=str(record.seq))
            if parsed is None:
                print(f"Warning: Could not parse header for sequence ID {seq_id}, skipping.")
                continue
            filtered_records.append(parsed)
        else:
            print(f"Warning: Sequence ID {seq_id} not found in FASTA file.")

    print(f"Writing {len(filtered_records)} sequences to {args.output_fasta}...")

    write_fasta_output(filtered_records, args.output_fasta)

    print(f"Writing sequence IDs to {args.output_ids}...")

    with open(args.output_ids, 'w') as f:
        json.dump(list(motif_hits), f, indent=2)


if __name__ == "__main__":
    main()