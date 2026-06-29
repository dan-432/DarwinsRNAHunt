#!/usr/bin/env python3
"""
Extract flanking CDS (coding sequences) around an RNA motif.

Given motif coordinates and a GFF annotation file, this script finds the
CDS features that flank the motif and extracts their protein sequences.
"""

import argparse
import sys
import tempfile
from pathlib import Path
import pandas
from DarwinsRNAHunt.genome_analysis import (
    load_genome_annotation,
    parse_fasta,
    get_cds_by_coordinates,
    extract_protein_id_from_gff,
    extract_cds_protein,
    write_fasta_output
)

COLUMNS=["seq_id", "start", "end", "bit_score"]

def main():
    parser = argparse.ArgumentParser(
        description="Extract flanking CDS around an RNA motif"
    )
    parser.add_argument(
        '--regions',
        required=True,
        type=str,
        help='Table of motif coordinates as tsv. "seq_id"\\t"start"\\t"end"\\t"bit_score"'
    )
    parser.add_argument(
        '--family',
        required=True,
        help='Rfam Family id, e.g. RF00050'
    )
    parser.add_argument(
        '--gff',
        required=True,
        help='GFF3 annotation file'
    )
    parser.add_argument(
        '--proteome',
        required=True,
        help='Protein FASTA file'
    )
    parser.add_argument(
        '--net-size',
        required=True,
        type=int,
        help='Nucleotides upand downstream of motif to grab cds'
    )
    parser.add_argument(
        '--assembly-accession',
        required=True,
        help='genome'
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Output protein FASTA file'
    )

    args = parser.parse_args()

    print(f"Parsing regions table from {args.regions}...", file=sys.stderr)

    motif_tbl = pandas.read_table(args.regions, header=None, names=COLUMNS)

    print(f"Loading proteome from {args.proteome}", file=sys.stderr)
    proteins = parse_fasta(args.proteome)

    print(f"Loading GFF annotation from {args.gff}", file=sys.stderr)
    with tempfile.NamedTemporaryFile(suffix='.db', delete=False) as tmp_db:
        db_path = tmp_db.name

    annotation_db = load_genome_annotation(args.gff, db_path)

    output_prot_seqs = []

    for _, motif in motif_tbl.iterrows():
        seq_id = motif["seq_id"]
        start = motif["start"]
        end = motif["end"]
        bit_score = motif["bit_score"]
        
        print(f"Processing motif {start}-{end} in {seq_id} (score: {bit_score})")
        
        surrounding_cds = get_cds_by_coordinates(annotation_db, seq_id, start, end, int(args.net_size))
        
        for cds in surrounding_cds:
            protein_id = extract_protein_id_from_gff(cds.attributes)
            if protein_id:
                seq = extract_cds_protein(protein_id, proteins)
                if seq:
                    output_prot_seqs.append({
                        "protein_id": f"{args.assembly_accession}|{protein_id}",
                        "type": "Protein",
                        "seqid": cds.seqid,
                        "start": cds.start,
                        "end": cds.end,
                        "strand": cds.strand,
                        "sequence": seq
                    })
                    print(f"Added flanking CDS: {protein_id}", file=sys.stderr)
                else:
                    print(f"Warning: CDS protein not found in proteome: {protein_id}", file=sys.stderr)
            else:
                print(f"Warning: Could not extract protein ID from CDS", file=sys.stderr)

    if output_prot_seqs:
        write_fasta_output(output_prot_seqs, args.output)
        print(f"Wrote {len(output_prot_seqs)} sequences to {args.output}", file=sys.stderr)
    else:
        print("Warning: No flanking CDS proteins to output, creating empty file", file=sys.stderr)
        Path(args.output).touch()

    # Cleanup temporary database
    try:
        Path(db_path).unlink()
    except:
        pass     


if __name__ == '__main__':
    main()

