#!/usr/bin/env python3
"""
Extract fixed flanking regions (upstream/downstream) around genes containing protein domains.
Simple approach: gene start - upstream_bp to gene end + downstream_bp
"""

import argparse
import os
import tempfile
import sys
from Bio import SeqIO, SearchIO
from DarwinsRNAHunt.genome_analysis import extract_flanking, parse_hmmer_hits, load_genome_annotation, parse_fasta, write_fasta_output, write_bed_output

# def parse_hmmer_hits(tblout_file, evalue_threshold=1e-5):
#     """Parses HMMER3 tblout file and returns a set of unique target IDs.
#     Args: tblout_file: path to HMMER3 tblout file
#           evalue_threshold: maximum e-value to consider a hit valid (default: 1e-5)
#     Returns: set of target IDs with hits below the e-value threshold"""
#     hit_ids = set()
#     # Use 'hmmer3-tab' for --tblout, 'hmmer3-domtab' for --domtblout
#     for query_result in SearchIO.parse(tblout_file, 'hmmer3-tab'):
#         for hit in query_result.hits:
#             if hit.evalue <= evalue_threshold:
#                 hit_ids.add(hit.id)
#     return hit_ids


# def get_cds_from_gff(annotation_db, hmmer_hit_ids):
#     """
#     CAN BE OPTIMISED - DON'T REALLY NEED TO CREATE A FULL DATABASE, JUST NEED TO ITERATE THROUGH CDS FEATURES AND CHECK FOR MATCHES
#     Matches HMMER hits to gene features in a GFF3 file.
#     Args: annotation_db: gffutils database object with genome annotation loaded
#           hmmer_hit_ids: set of protein IDs from HMMER hits
#     Returns: list of gffutils.Feature objects for matched genes"""
    
#     matched_genes = []
        
#     # Iterate through hits and find corresponding genes in GFF
#     for hit_id in hmmer_hit_ids:
#         try:
#             # look for cds feature with matching protein_id
#             feature = annotation_db["cds-" + hit_id]
#             # If hit is not gene itself, get parent gene feature
#             if feature.featuretype != 'CDS':
#                 print(f"Warning: Feature {feature.id} is not a CDS, skipping", file=sys.stderr)
#                 continue

#             matched_genes.append(feature)
#         except gffutils.exceptions.FeatureNotFoundError:
#             continue
    
#     return matched_genes


# def extract_segment(seq, start, end, strand):
#     """Extracts a segment from a sequence and reverse complements if on minus strand.
#     Args: seq: Bio.Seq object
#           start: 0-based start position
#           end: 0-based end position
#           strand: '+' or '-'
#     Returns: extracted segment as Bio.Seq object"""
#     segment = seq[start:end]
#     return segment.reverse_complement() if strand == '-' else segment

# def extract_gene_regions(genome_file, gene_coords, upstream_bp, downstream_bp):
#     """Extract 5' UTR, CDS, and 3' UTR regions for each gene.
#     Args: genome_file: path to genome FASTA
#           gene_coords: list of gffutils.Feature objects for matched genes
#           upstream_bp: dictionary {protein_id: number of bases to extract upstream of gene start}
#           downstream_bp: dictionary {protein_id: number of bases to extract downstream of gene end}
#     Returns: list of dicts with sequence info for utr5, cds, and utr3"""
#     genome_seqs = SeqIO.to_dict(SeqIO.parse(genome_file, 'fasta'))
#     extracted_regions = []

#     for feature in gene_coords:
#         seqid = feature.chrom
#         strand = feature.strand
#         gene_start = int(feature.start)
#         gene_end = int(feature.end)

#         if seqid not in genome_seqs:
#             print(f"Warning: Sequence {seqid} not found in genome", file=sys.stderr)
#             continue

#         seq = genome_seqs[seqid].seq
#         py_start = gene_start - 1
#         py_end = gene_end

#         if strand == '+':
#             utr5_start = max(0, py_start - upstream_bp[feature.id])
#             utr5_end = py_start
#             utr3_start = py_end
#             utr3_end = min(len(seq), py_end + downstream_bp[feature.id])
#         else:
#             utr5_start = py_end
#             utr5_end = min(len(seq), py_end + upstream_bp[feature.id])
#             utr3_start = max(0, py_start - downstream_bp[feature.id])
#             utr3_end = py_start

#         utr5_seq = extract_segment(seq, utr5_start, utr5_end, strand)
#         cds_seq = extract_segment(seq, py_start, py_end, strand)
#         utr3_seq = extract_segment(seq, utr3_start, utr3_end, strand)

#         extracted_regions.extend([
#             {
#                 'protein_id': feature.id,
#                 'type': "utr5",
#                 'seqid': seqid,
#                 'start': utr5_start + 1,
#                 'end': utr5_end,
#                 'strand': strand,
#                 'gene_start': gene_start,
#                 'gene_end': gene_end,
#                 'sequence': utr5_seq
#             },
#             {
#                 'protein_id': feature.id,
#                 'type': "cds",
#                 'seqid': seqid,
#                 'start': py_start + 1,
#                 'end': py_end,
#                 'strand': strand,
#                 'gene_start': gene_start,
#                 'gene_end': gene_end,
#                 'sequence': cds_seq
#             },
#             {
#                 'protein_id': feature.id,
#                 'type': "utr3",
#                 'seqid': seqid,
#                 'start': utr3_start + 1,
#                 'end': utr3_end,
#                 'strand': strand,
#                 'gene_start': gene_start,
#                 'gene_end': gene_end,
#                 'sequence': utr3_seq
#             }
#         ])

#     return extracted_regions

# def get_next_gene_distance(annotation_db, genes, upstream=True):
#     """Calculate distance to next gene in specified direction for each gene.
#     Args: annotation_db: gffutils database object
#           genes: list of gffutils.Feature objects for matched genes
#           upstream: True for upstream, False for downstream
#     Returns: dict {protein_id: distance to next gene in specified direction}"""
    
#     distances = {}
    
#     for feature in genes:
#         direcion = "<="
#         if upstream:
#             if feature.strand == '+':
#                 direcion = "<="
#             else:
#                 direcion = ">="
#         else:
#             if feature.strand == '+':
#                 direcion = ">="
#             else:
#                 direcion = "<="

#         # Find nearest neighbor in specified direction
#         query = f"""
#                 SELECT * FROM features
#                 WHERE seqid = '{feature.seqid}'
#                 AND featuretype = 'CDS'
#                 AND start {direcion} {feature.start if upstream else feature.end}
#                 AND id != '{feature.id}'
#                 ORDER BY ABS(start - {feature.start if upstream else feature.end})
#                 LIMIT 1
#                 """
#         neighbor = list(annotation_db.execute(query))
#         if not neighbor:
#             distances[feature.id] = float('inf')
#             continue

#         # note this give sqlite3.row object, not gffutils.Feature, so we need to convert it back to a Feature
#         neighbor = neighbor[0]
#         # Calculate gap
#         # Return absolute distance between gene and neighbor, ensuring it's non-negative and being strand-aware
#         if feature.strand == '+':
#             if upstream:
#                 distance = max(0, feature.start - neighbor["end"])
#             else:
#                 distance = max(0, neighbor["start"] - feature.end)
#         else:
#             if upstream:
#                 distance = max(0, neighbor["start"] - feature.end)
#             else:
#                 distance = max(0, feature.start - neighbor["end"])

#         print("Distance for feature {} to neighbor {}: {}".format(feature.id, neighbor["id"], distance), file=sys.stderr)

#         distances[feature.id] = distance

#     return distances

        


# def write_fasta_output(regions, utr5_path, cds_path, utr3_path):
#     """Write extracted regions to separate FASTA files.
#     Args: regions: list of dicts with sequence info
#           utr5_path: path to output FASTA for 5' UTRs
#           cds_path: path to output FASTA for CDS
#           utr3_path: path to output FASTA for 3' UTRs"""
#     with open(utr5_path, "w") as f5, \
#          open(cds_path, "w") as fcds, \
#          open(utr3_path, "w") as f3:
#         for region in regions:
#             header = (
#                 f">{region['protein_id']}|{region['type']}|"
#                 f"{region['seqid']}:{region['start']}-{region['end']}({region['strand']})"
#             )
#             out = {
#                 "utr5": f5,
#                 "cds": fcds,
#                 "utr3": f3
#             }[region["type"]]
#             out.write(header + "\n")
#             out.write(str(region["sequence"]) + "\n")

# def write_bed_output(regions, output_file):
#     """Write extracted regions to BED format.
#     Args: regions: list of dicts with sequence info
#           output_file: path to output BED file"""
#     with open(output_file, 'w') as f:
#         f.write("# chrom\tstart\tend\tname\tscore\tstrand\ttype\tgene_start\tgene_end\n")
#         for region in regions:
#             bed_start = region['start'] - 1
#             bed_end = region['end']
#             name = f"{region['protein_id']}|{region['type']}"
#             f.write(
#                 f"{region['seqid']}\t{bed_start}\t{bed_end}\t{name}\t0\t"
#                 f"{region['strand']}\t{region['type']}\t{region.get('gene_start', '')}\t{region.get('gene_end', '')}\n"
#             )
    
# def load_genome_annotation(gff_file, db_path):
#     """Load genome annotation into a database.
#     Args: gff_file: path to GFF file
#           db_path: path to the database file
#     Returns: gffutils.FeatureDB object with genome annotation loaded"""
    
#     db = gffutils.create_db(gff_file, dbfn=db_path, force=True, 
#                             keep_order=True, merge_strategy='merge', 
#                             sort_attribute_values=True)
#     return db



# def pmain():
#     print("EXTRACTING FLANKING SEQUENCES")

#     parser = argparse.ArgumentParser(
#         description='Extract fixed flanking regions around domain-containing genes'
#     )
#     parser.add_argument('--input-hits', required=True,
#                        help='HMMER hits file (tsv format from --tblout)')
#     parser.add_argument('--genome', required=True,
#                        help='Genome FASTA file')
#     parser.add_argument('--annotation', required=True,
#                        help='GFF annotation file')
#     parser.add_argument('--output-upstream', required=True,
#                         help='Output FASTA file for upstream flanking sequences')
#     parser.add_argument('--output-cds', required=True,
#                         help='Output FASTA file for CDS sequences')
#     parser.add_argument('--output-downstream', required=True,
#                         help='Output FASTA file for downstream flanking sequences')
#     parser.add_argument('--output-coords', required=True,
#                         help='Output BED file')
#     parser.add_argument('--upstream', type=int, default=500,
#                        help='Bases upstream of gene start, -1 for full intergenic sequence (default: 500)')
#     parser.add_argument('--downstream', type=int, default=500,
#                        help='Bases downstream of gene end, -1 for full intergenic sequence (default: 500)')
#     parser.add_argument('--evalue', type=float, default=1e-5,
#                        help='E-value threshold (default: 1e-5)')
#     # example usage:
#     # python workflow/src/extract_flanking.py --input-hits results/domain_hits/GCA_000013605.1_protdomain_hits.tsv --genome resources/genomes/GCA_000013605.1/genomic.fna --annotation resources/genomes/GCA_000013605.1/geneannotation.gff --output-upstream results/utr5.fasta --output-cds results/cds.fasta --output-downstream results/utr3.fasta --output-coords results/coords.bed --upstream -1 --downstream -1 --evalue 1e-5

#     args = parser.parse_args()
    
#     print("Parsing HMMER hits...", file=sys.stderr)
#     protein_hits = parse_hmmer_hits(args.input_hits, args.evalue)
#     print(f"Found {len(protein_hits)} proteins with domain hits", file=sys.stderr)

#     # Use a unique temporary database file to avoid conflicts in parallel runs
#     with tempfile.NamedTemporaryFile(suffix='.db', delete=False, dir=".temp") as tmp_db:
#         db_path = tmp_db.name
#     annotation_db = load_genome_annotation(args.annotation, db_path)
    
#     print(f"Extracting CDS coordinates from {args.annotation}...", file=sys.stderr)
#     gene_coords = get_cds_from_gff(annotation_db, protein_hits)
#     print(f"Found coordinates for {len(gene_coords)} CDSs", file=sys.stderr)
    
    
#     # make dictionaries of flanking distances if they are set to a specific value
#     up_flanking = {feature.id: args.upstream for feature in gene_coords}
#     down_flanking = {feature.id: args.downstream for feature in gene_coords}
    
#     if args.upstream == -1 or args.downstream == -1:
#         print("Warning: Using full intergenic regions for flanking sequences", file=sys.stderr)

#     if args.upstream == -1:
#         print("Calculating full intergenic upstream regions...", file=sys.stderr)

#         up_flanking = get_next_gene_distance(annotation_db, gene_coords, upstream=True)
#         print(f"Calculated upstream flanking distances: {up_flanking}", file=sys.stderr)
#     if args.downstream == -1:
#         print("Calculating full intergenic downstream regions...", file=sys.stderr)

#         down_flanking = get_next_gene_distance(annotation_db, gene_coords, upstream=False)
#         print(f"Calculated downstream flanking distances: {down_flanking}", file=sys.stderr)

#     genome_seqs = SeqIO.to_dict(SeqIO.parse(args.genome, 'fasta'))

#     regions = extract_gene_regions(
#         genome_seqs,
#         gene_coords,
#         up_flanking,
#         down_flanking
#     )

#     write_fasta_output(
#         regions,
#         args.output_upstream,
#         args.output_cds,
#         args.output_downstream
#     )

#     write_bed_output(
#         regions,
#         args.output_coords
#     )
    
#     print(f"Wrote sequences to {args.output_upstream}, {args.output_cds}, {args.output_downstream}", file=sys.stderr)
#     print(f"Wrote coordinates to {args.output_coords}", file=sys.stderr)

#     print("Removing temporary database...", file=sys.stderr)

#     # Ensure the temporary database is removed
#     try:
#         os.unlink(db_path)
#     except OSError:
#         pass

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
    # example usage:
    # python workflow/src/extract_flanking.py --input-hits results/domain_hits/GCA_000013605.1_protdomain_hits.tsv --genome resources/genomes/GCA_000013605.1/genomic.fna --annotation resources/genomes/GCA_000013605.1/geneannotation.gff --output-upstream results/utr5.fasta --output-cds results/cds.fasta --output-downstream results/utr3.fasta --output-coords results/coords.bed --upstream -1 --downstream -1 --evalue 1e-5

    args = parser.parse_args()
    
    print("Parsing HMMER hits...", file=sys.stderr)
    protein_hits = parse_hmmer_hits(args.input_hits, args.evalue)
    print(f"Found {len(protein_hits)} proteins with domain hits", file=sys.stderr)

    # Use a unique temporary database file to avoid conflicts in parallel runs
    with tempfile.NamedTemporaryFile(suffix='.db', delete=False, dir=".temp") as tmp_db:
        db_path = tmp_db.name
    annotation_db = load_genome_annotation(args.annotation, db_path)

    genome_seq = parse_fasta(args.genome)
    
    regions = extract_flanking(protein_hits, annotation_db, args.upstream, args.downstream, genome_seq)

    write_fasta_output(
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