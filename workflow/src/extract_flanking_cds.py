"""
Extract flanking CDS (coding sequences) around an RNA motif.

Given motif coordinates and a GFF annotation file, this script finds the
CDS features that flank the motif and extracts their protein sequences.
It also emits a position table recording where each extracted CDS sits
relative to the motif (position, rank, distance, strand), so that this
information can be joined back in after clustering.
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
    extract_cds_protein,
    write_fasta_output,
)
from DarwinsRNAHunt.sequence import ProteinSequence
from DarwinsRNAHunt.motif import Motif, POSITION_COLUMNS

COLUMNS = ["seq_id", "start", "end", "bit_score"]


def main():
    parser = argparse.ArgumentParser(
        description="Extract flanking CDS around an RNA motif"
    )
    parser.add_argument('--regions', required=True, type=str,
        help='Table of motif coordinates as tsv. "seq_id"\\t"start"\\t"end"\\t"bit_score"')
    parser.add_argument('--family', required=True, help='Rfam Family id, e.g. RF00050')
    parser.add_argument('--gff', required=True, help='GFF3 annotation file')
    parser.add_argument('--proteome', required=True, help='Protein FASTA file')
    parser.add_argument('--net-size', required=True, type=int,
        help='Nucleotides upand downstream of motif to grab cds')
    parser.add_argument('--assembly-accession', required=True, help='genome')
    parser.add_argument('--output', required=True, help='Output protein FASTA file')
    parser.add_argument('--position-table', required=True,
        help='Table of motif positions relative to the motif to include in output')

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
    position_rows = []

    for _, motif_row in motif_tbl.iterrows():
        seq_id = motif_row["seq_id"]
        bit_score = motif_row["bit_score"]

        # motifs.tsv has no explicit strand column -- strand is encoded
        # implicitly by coordinate order (cmsearch reports "seq from" >
        # "seq to" for minus-strand hits). Motif.from_hit derives strand
        # from raw start/end order and normalizes coordinates itself.
        raw_start, raw_end = motif_row["start"], motif_row["end"]
        start, end = min(raw_start, raw_end), max(raw_start, raw_end)

        motif = Motif.from_hit(
            id=f"{seq_id}:{start}-{end}",
            raw_start=raw_start,
            raw_end=raw_end,
            sequence="",   # sequence not available from --regions; unused for positions
            seq_id=seq_id,
        )

        print(f"Processing motif {start}-{end} in {seq_id} "
              f"(strand: {motif.strand.value}, score: {bit_score})")

        surrounding_cds = get_cds_by_coordinates(annotation_db, seq_id, start, end, int(args.net_size))

        # Strand-aware classification + within-subset ranking, done in one pass.
        cds_positions = motif.get_cds_positions_relative_to(surrounding_cds)

        if not motif.fully_intergenic:
            print(f"Note: motif {start}-{end} in {seq_id} overlaps at least "
                  f"one flanking CDS", file=sys.stderr)

        for cds_pos in cds_positions:
            cds = cds_pos.cds  # CDSSequence, seq=None (no genome FASTA in get_cds_by_coordinates)

            protein_id = cds.get_protein_id()
            if not protein_id:
                print(f"Warning: Could not extract protein ID from CDS", file=sys.stderr)
                continue

            seq = extract_cds_protein(protein_id, proteins)
            if not seq:
                print(f"Warning: CDS protein not found in proteome: {protein_id}", file=sys.stderr)
                continue

            protein = ProteinSequence(
                id=protein_id,
                seq=seq,
                seqid=cds.seqid,
                start=cds.start,
                end=cds.end,
                strand=cds.strand,
                assembly_accession=args.assembly_accession,
            )
            output_prot_seqs.append(protein)

            # header doubles as the join key downstream (kept from the
            # original script's approach -- gross, but it works)
            cds_pos.protein = protein

            print(f"Added flanking CDS: {protein_id} "
                  f"({cds_pos.position.value}, rank {cds_pos.rank}, "
                  f"{cds_pos.distance} bp)", file=sys.stderr)

        position_rows.extend(motif.to_rows())

    if output_prot_seqs:
        write_fasta_output(output_prot_seqs, args.output)
        print(f"Wrote {len(output_prot_seqs)} sequences to {args.output}", file=sys.stderr)
    else:
        print("Warning: No flanking CDS proteins to output, creating empty file", file=sys.stderr)
        Path(args.output).touch()

    output_position_tbl = pandas.DataFrame(position_rows, columns=POSITION_COLUMNS)
    output_position_tbl.to_csv(args.position_table, sep="\t", index=False)
    print(f"Wrote {len(output_position_tbl)} rows to {args.position_table}", file=sys.stderr)

    try:
        Path(db_path).unlink()
    except:
        pass


if __name__ == '__main__':
    main()