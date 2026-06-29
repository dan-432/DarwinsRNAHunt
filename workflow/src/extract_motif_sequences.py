import argparse
import json
from DarwinsRNAHunt.genome_analysis import parse_tblout, parse_fasta, write_fasta_output

def main():
    # example usage:
    # python workflow/src/extract_motif_sequences.py --fasta-file results/02_sequence_selection/03_combined/{target_domain}_domain_flanking_upstream.fasta --motif-hits-tbl results/03_motif_discovery/{target_domain}/02_motif_hits/{motif_id}_hits.tblout --output-fasta results/03_motif_discovery/{target_domain}/03_homolog_search/expanded_seqs/{motif_id}_expanded.fasta --output-ids results/03_motif_discovery/{target_domain}/03_homolog_search/expanded_seqs/{motif_id}_ids.json
    parser = argparse.ArgumentParser(description="Extract motif sequences from a FASTA file.")
    parser.add_argument("--fasta-file", help="Path to the input FASTA file.")
    parser.add_argument("--motif-hits-tbl", help="Path to the motif hits table.")
    parser.add_argument("--output-fasta", help="Path to the output FASTA file.")
    parser.add_argument("--output-ids", help="Path to the output file for sequence IDs.")
    
    args = parser.parse_args()
    
    seqs = parse_fasta(args.fasta_file)
    motif_hits = parse_tblout(args.motif_hits_tbl)

    filtered_records = [seqs[id] for id in motif_hits if id in seqs]

    write_fasta_output(filtered_records, args.output_fasta)

    json.dump(list(motif_hits), open(args.output_ids, 'w'), indent=2)