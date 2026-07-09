#!/usr/bin/env python3
"""
Functions for extracting sequences and basic gff and fasta parsing
"""

import sys
from Bio import SeqIO, SearchIO
import gffutils

def parse_hmmer_hits(tblout_file, evalue_threshold=1e-5):
    """Parses HMMER3 tblout file and returns a set of unique target IDs.
    Args: tblout_file: path to HMMER3 tblout file
          evalue_threshold: maximum e-value to consider a hit valid (default: 1e-5)
    Returns: set of target IDs with hits below the e-value threshold"""
    hit_ids = set()
    # Use 'hmmer3-tab' for --tblout, 'hmmer3-domtab' for --domtblout
    for query_result in SearchIO.parse(tblout_file, 'hmmer3-tab'):
        for hit in query_result.hits:
            if hit.evalue <= evalue_threshold:
                hit_ids.add(hit.id)
    return hit_ids

def parse_fasta(fasta_file):
    """Parses fasta file to dictionary of Seq IO records, deduplicates by sequence ID."""
    seq_dict = {}
    for record in SeqIO.parse(fasta_file, 'fasta'):
        seq_dict[record.id] = record
    return seq_dict

def get_cds_from_gff(annotation_db, hmmer_hit_ids):
    """
    CAN BE OPTIMISED - DON'T REALLY NEED TO CREATE A FULL DATABASE, JUST NEED TO ITERATE THROUGH CDS FEATURES AND CHECK FOR MATCHES
    Matches HMMER hits to gene features in a GFF3 file.
    Args: annotation_db: gffutils database object with genome annotation loaded
          hmmer_hit_ids: set of protein IDs from HMMER hits
    Returns: list of gffutils.Feature objects for matched genes"""
    
    matched_genes = []
        
    # Iterate through hits and find corresponding genes in GFF
    for hit_id in hmmer_hit_ids:
        try:
            # look for cds feature with matching protein_id
            feature = annotation_db["cds-" + hit_id]
            # If hit is not gene itself, get parent gene feature
            if feature.featuretype != 'CDS':
                print(f"Warning: Feature {feature.id} is not a CDS, skipping", file=sys.stderr)
                continue

            matched_genes.append(feature)
        except gffutils.exceptions.FeatureNotFoundError:
            continue
    
    return matched_genes


def extract_segment(seq, start, end, strand):
    """Extracts a segment from a sequence and reverse complements if on minus strand.
    Args: seq: Bio.Seq object
          start: 0-based start position
          end: 0-based end position
          strand: '+' or '-'
    Returns: extracted segment as Bio.Seq object"""
    segment = seq[start:end]
    return segment.reverse_complement() if strand == '-' else segment

def extract_gene_regions(genome_seqs, gene_coords, upstream_bp, downstream_bp):
    """Extract 5' UTR, CDS, and 3' UTR regions for each gene.
    Args: genome_seqs: FASTA parsed to dictionary
          gene_coords: list of gffutils.Feature objects for matched genes
          upstream_bp: dictionary {protein_id: number of bases to extract upstream of gene start}
          downstream_bp: dictionary {protein_id: number of bases to extract downstream of gene end}
    Returns: list of dicts with sequence info for utr5, cds, and utr3"""
    extracted_regions = []

    for feature in gene_coords:
        seqid = feature.chrom
        strand = feature.strand
        gene_start = int(feature.start)
        gene_end = int(feature.end)

        if seqid not in genome_seqs:
            print(f"Warning: Sequence {seqid} not found in genome", file=sys.stderr)
            continue

        seq = genome_seqs[seqid].seq
        py_start = gene_start - 1
        py_end = gene_end

        if strand == '+':
            utr5_start = max(0, py_start - upstream_bp[feature.id])
            utr5_end = py_start
            utr3_start = py_end
            utr3_end = min(len(seq), py_end + downstream_bp[feature.id])
        else:
            utr5_start = py_end
            utr5_end = min(len(seq), py_end + upstream_bp[feature.id])
            utr3_start = max(0, py_start - downstream_bp[feature.id])
            utr3_end = py_start

        utr5_seq = extract_segment(seq, utr5_start, utr5_end, strand)
        cds_seq = extract_segment(seq, py_start, py_end, strand)
        utr3_seq = extract_segment(seq, utr3_start, utr3_end, strand)

        extracted_regions.extend([
            {
                'protein_id': feature.id,
                'type': "utr5",
                'seqid': seqid,
                'start': utr5_start + 1,
                'end': utr5_end,
                'strand': strand,
                'gene_start': gene_start,
                'gene_end': gene_end,
                'sequence': utr5_seq
            },
            {
                'protein_id': feature.id,
                'type': "cds",
                'seqid': seqid,
                'start': py_start + 1,
                'end': py_end,
                'strand': strand,
                'gene_start': gene_start,
                'gene_end': gene_end,
                'sequence': cds_seq
            },
            {
                'protein_id': feature.id,
                'type': "utr3",
                'seqid': seqid,
                'start': utr3_start + 1,
                'end': utr3_end,
                'strand': strand,
                'gene_start': gene_start,
                'gene_end': gene_end,
                'sequence': utr3_seq
            }
        ])

    return extracted_regions

def get_next_gene_distance(annotation_db, genes, upstream=True):
    """Calculate distance to next gene in specified direction for each gene.
    Args: annotation_db: gffutils database object
          genes: list of gffutils.Feature objects for matched genes
          upstream: True for upstream, False for downstream
    Returns: dict {protein_id: distance to next gene in specified direction}"""
    
    distances = {}
    
    for feature in genes:
        direcion = "<="
        if upstream:
            if feature.strand == '+':
                direcion = "<="
            else:
                direcion = ">="
        else:
            if feature.strand == '+':
                direcion = ">="
            else:
                direcion = "<="

        # Find nearest neighbor in specified direction
        query = f"""
                SELECT * FROM features
                WHERE seqid = '{feature.seqid}'
                AND featuretype = 'CDS'
                AND start {direcion} {feature.start if upstream else feature.end}
                AND id != '{feature.id}'
                ORDER BY ABS(start - {feature.start if upstream else feature.end})
                LIMIT 1
                """
        neighbor = list(annotation_db.execute(query))
        if not neighbor:
            distances[feature.id] = float('inf')
            continue

        # note this give sqlite3.row object, not gffutils.Feature, so we need to convert it back to a Feature
        neighbor = neighbor[0]
        # Calculate gap
        # Return absolute distance between gene and neighbor, ensuring it's non-negative and being strand-aware
        if feature.strand == '+':
            if upstream:
                distance = max(0, feature.start - neighbor["end"])
            else:
                distance = max(0, neighbor["start"] - feature.end)
        else:
            if upstream:
                distance = max(0, neighbor["start"] - feature.end)
            else:
                distance = max(0, feature.start - neighbor["end"])

        print("Distance for feature {} to neighbor {}: {}".format(feature.id, neighbor["id"], distance), file=sys.stderr)

        distances[feature.id] = distance

    return distances


def write_fasta_output(regions, output_file):
    """Write extracted regions to a single FASTA file.
    Args: regions: list of dicts with sequence info
          output_file: path to output FASTA file"""
    with open(output_file, "w") as f:
        for region in regions:
            header = (
                f">{region['protein_id']}|{region['type']}|"
                f"{region['seqid']}:{region['start']}-{region['end']}({region['strand']})"
            )
            f.write(header + "\n")
            f.write(str(region["sequence"]) + "\n")        


def write_fasta_outputs(regions, utr5_path, cds_path, utr3_path):
    """Write extracted regions to separate FASTA files.
    Args: regions: list of dicts with sequence info
          utr5_path: path to output FASTA for 5' UTRs
          cds_path: path to output FASTA for CDS
          utr3_path: path to output FASTA for 3' UTRs"""
    write_fasta_output([r for r in regions if r['type'] == 'utr5'], utr5_path)
    write_fasta_output([r for r in regions if r['type'] == 'cds'], cds_path)
    write_fasta_output([r for r in regions if r['type'] == 'utr3'], utr3_path)


def write_bed_output(regions, output_file):
    """Write extracted regions to BED format.
    Args: regions: list of dicts with sequence info
          output_file: path to output BED file"""
    with open(output_file, 'w') as f:
        f.write("# chrom\tstart\tend\tname\tscore\tstrand\ttype\tgene_start\tgene_end\n")
        for region in regions:
            bed_start = region['start'] - 1
            bed_end = region['end']
            name = f"{region['protein_id']}|{region['type']}"
            f.write(
                f"{region['seqid']}\t{bed_start}\t{bed_end}\t{name}\t0\t"
                f"{region['strand']}\t{region['type']}\t{region.get('gene_start', '')}\t{region.get('gene_end', '')}\n"
            )
    
def load_genome_annotation(gff_file, db_path):
    """Load genome annotation into a database.
    Args: gff_file: path to GFF file
          db_path: path to the database file
    Returns: gffutils.FeatureDB object with genome annotation loaded"""
    
    db = gffutils.create_db(gff_file, dbfn=db_path, force=True, 
                            keep_order=True, merge_strategy='merge', 
                            sort_attribute_values=True)
    return db



def extract_flanking(protein_hits, annotation_db, upstream: int, downstream: int, genome_seqs: dict):
    """
    Extracts flanking sequences for genes corresponding to HMMER hits.
    Args: protein_hits: set of protein IDs from HMMER hits
      annotation_db: gffutils database object with genome annotation loaded
      upstream: number of bases to extract upstream of gene start (set to -1 to use full intergenic region)
      downstream: number of bases to extract downstream of gene end (set to -1 to use full intergenic region)
      genome_seqs: dictionary of genome sequences parsed from FASTA file
    Returns: list of dicts with sequence info for utr5, cds, and utr3
    """

    gene_coords = get_cds_from_gff(annotation_db, protein_hits)

    print(f"Found coordinates for {len(gene_coords)} CDSs", file=sys.stderr)
    
    # make dictionaries of flanking distances if they are set to a specific value
    up_flanking = {feature.id: upstream for feature in gene_coords}
    down_flanking = {feature.id: downstream for feature in gene_coords}
    
    if upstream == -1 or downstream == -1:
        print("Warning: Using full intergenic regions for flanking sequences", file=sys.stderr)

    if upstream == -1:
        print("Calculating full intergenic upstream regions...", file=sys.stderr)

        up_flanking = get_next_gene_distance(annotation_db, gene_coords, upstream=True)
        print(f"Calculated upstream flanking distances: {up_flanking}", file=sys.stderr)
    if downstream == -1:
        print("Calculating full intergenic downstream regions...", file=sys.stderr)

        down_flanking = get_next_gene_distance(annotation_db, gene_coords, upstream=False)
        print(f"Calculated downstream flanking distances: {down_flanking}", file=sys.stderr)

    regions = extract_gene_regions(
        genome_seqs,
        gene_coords,
        up_flanking,
        down_flanking
    )

    return regions

def parse_tblout(tblout_file):
    """Parse motif hits from cmfinder tblout file.
    Args: tblout_file: path to cmfinder tblout file
    Returns: set of sequence IDs with motif hits"""
    motif_hits = set()
    with open(tblout_file, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                fields = line.strip().split()
                if len(fields) >= 2:
                    motif_hits.add(fields[0])
    return motif_hits


def get_cds_by_coordinates(annotation_db, seqid, motif_start, motif_end, net_size):
    """Find CDS features flanking motif coordinates.
    Args: annotation_db: gffutils database object with genome annotation loaded
          seqid: sequence/chromosome ID
          motif_start: motif start coordinate (1-based)
          motif_end: motif end coordinate (1-based)
          net_size: how far up and downstream to extract features
    Returns: set of neighbouring gffutils.Feature objects"""

    lower_bound = max(1, motif_start - net_size)
    upper_bound = motif_end + net_size
    
    # # Query for CDS features on this sequence
    # query = f"""
    #     SELECT * FROM features
    #     WHERE seqid = '{seqid}'
    #     AND featuretype = 'CDS'
    #     AND start <= {upper_bound}
    #     AND end >= {lower_bound}
    #     ORDER BY start
    # """

    # # collect anything in the bounds of our search net
    
    # cds_features = set(annotation_db.execute(query))
    
    # return cds_features

    # Query for CDS features on this sequence using gffutils region()
    cds_features = list(annotation_db.region(
        seqid=seqid,
        start=lower_bound,
        end=upper_bound,
        featuretype='CDS',
        completely_within=False
    ))
    
    return cds_features


def extract_cds_protein(protein_id, proteome):
    """Extract protein sequence by ID from proteome dictionary.
    Args: protein_id: protein identifier
          proteome: dictionary of protein sequences (from parse_fasta)
    Returns: protein sequence string or None if not found"""
    
    # Try exact match first
    if protein_id in proteome:
        return str(proteome[protein_id].seq)
    
    # Try partial match (in case IDs are truncated or have prefixes)
    for seq_id in proteome:
        if protein_id in seq_id or seq_id in protein_id:
            return str(proteome[seq_id].seq)
    
    return None


def extract_protein_id_from_gff(attributes):
    """Extract protein ID from GFF attributes dict.
    Args: attributes: dictionary of GFF attributes
    Returns: protein ID string or None"""
    
    # Try common attribute keys for protein IDs
    for key in ['protein_id', 'Name', 'ID', 'Dbxref']:
        if key in attributes:
            value = attributes[key][0]
            if key == 'Dbxref':
                # Handle format like "protein_id:XP_001234567"
                for item in value.split(','):
                    if 'protein_id:' in item or 'RefSeq:' in item:
                        return item.split(':')[1]
            else:
                return value
    
    return None