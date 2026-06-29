#!/usr/bin/env python3
"""
Extract fixed flanking regions (upstream/downstream) around genes containing protein domains.
Simple approach: gene start - upstream_bp to gene end + downstream_bp
"""

import argparse
import os
import tempfile
import sys
from DarwinsRNAHunt.genome_analysis import extract_flanking, parse_hmmer_hits, load_genome_annotation, parse_fasta, write_fasta_outputs, write_bed_output

def main():
    print("EXTRACTING FLANKING SEQUENCES")

    parser = argparse.ArgumentParser(
        description='Extract fixed flanking regions around domain-containing genes'
    )
    parser.add_argument('--input-hits', required=True,
                       help='HMMER hits file (tsv format from --tblout)')
    parser.add_argument('--genome', required=True,
                       help='Genome FASTA file')
    parser.add_argument('--annotation', required=True,
                       help='GFF annotation file')
    parser.add_argument('--output-upstream', required=True,
                        help='Output FASTA file for upstream flanking sequences')
    parser.add_argument('--output-cds', required=True,
                        help='Output FASTA file for CDS sequences')
    parser.add_argument('--output-downstream', required=True,
                        help='Output FASTA file for downstream flanking sequences')
    parser.add_argument('--output-coords', required=True,
                        help='Output BED file')
    parser.add_argument('--upstream', type=int, default=500,
                       help='Bases upstream of gene start, -1 for full intergenic sequence (default: 500)')
    parser.add_argument('--downstream', type=int, default=500,
                       help='Bases downstream of gene end, -1 for full intergenic sequence (default: 500)')
    parser.add_argument('--evalue', type=float, default=1e-5,
                       help='E-value threshold (default: 1e-5)')
    parser.add_argument('--temp-dir', default='.temp',
                       help='Directory for temporary files (default: .temp)')
    # example usage:
    # python workflow/src/extract_flanking.py --input-hits results/domain_hits/GCA_000013605.1_protdomain_hits.tsv --genome resources/genomes/GCA_000013605.1/genomic.fna --annotation resources/genomes/GCA_000013605.1/geneannotation.gff --output-upstream results/utr5.fasta --output-cds results/cds.fasta --output-downstream results/utr3.fasta --output-coords results/coords.bed --upstream -1 --downstream -1 --evalue 1e-5

    args = parser.parse_args()
    
    print("Parsing HMMER hits...", file=sys.stderr)
    protein_hits = parse_hmmer_hits(args.input_hits, args.evalue)
    print(f"Found {len(protein_hits)} proteins with domain hits", file=sys.stderr)

    # Use a unique temporary database file to avoid conflicts in parallel runs
    with tempfile.NamedTemporaryFile(suffix='.db', delete=False, dir=args.temp_dir) as tmp_db: # MAKE TEMP DIR ARGUMENT
        db_path = tmp_db.name
    annotation_db = load_genome_annotation(args.annotation, db_path)

    genome_seq = parse_fasta(args.genome)
    
    regions = extract_flanking(protein_hits, annotation_db, args.upstream, args.downstream, genome_seq)

    write_fasta_outputs(
        regions,
        args.output_upstream,
        args.output_cds,
        args.output_downstream
    )

    write_bed_output(
        regions,
        args.output_coords
    )
    
    print(f"Wrote sequences to {args.output_upstream}, {args.output_cds}, {args.output_downstream}", file=sys.stderr)
    print(f"Wrote coordinates to {args.output_coords}", file=sys.stderr)

    print("Removing temporary database...", file=sys.stderr)

    # Ensure the temporary database is removed
    try:
        os.unlink(db_path)
    except OSError:
        pass


if __name__ == '__main__':
    main()