"""
Functions for extracting sequences and basic gff and fasta parsing.

Sequence data is now modeled via the Sequence/NucleotideSequence/
ProteinSequence class hierarchy in DarwinsRNAHunt.sequence, rather than
hand-built dicts with a 'type' string key. Header encoding/decoding and
FASTA/BED writing are now methods on those classes (to_header, to_fasta,
to_bed_line) instead of free functions here.
"""

import sys
from Bio import SeqIO, SearchIO
import gffutils

from DarwinsRNAHunt.sequence import (
    Strand,
    NucleotideSequence,
    CDSSequence,
    UTR5Sequence,
    UTR3Sequence,
    ProteinSequence,
)


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

    for hit_id in hmmer_hit_ids:
        try:
            feature = annotation_db["cds-" + hit_id]
            if feature.featuretype != 'CDS':
                print(f"Warning: Feature {feature.id} is not a CDS, skipping", file=sys.stderr)
                continue

            matched_genes.append(feature)
        except gffutils.exceptions.FeatureNotFoundError:
            continue

    return matched_genes


def extract_gene_regions(genome_seqs, gene_coords, upstream_bp, downstream_bp):
    """Extract 5' UTR, CDS, and 3' UTR regions for each gene.
    Args: genome_seqs: FASTA parsed to dictionary (id -> SeqRecord)
          gene_coords: list of gffutils.Feature objects for matched genes
          upstream_bp: dictionary {protein_id: number of bases to extract upstream of gene start}
          downstream_bp: dictionary {protein_id: number of bases to extract downstream of gene end}
    Returns: list of NucleotideSequence subclass instances (UTR5Sequence, CDSSequence, UTR3Sequence)"""
    extracted_regions = []

    for feature in gene_coords:
        seqid = feature.chrom
        strand = Strand(feature.strand)
        gene_start = int(feature.start)
        gene_end = int(feature.end)

        if seqid not in genome_seqs:
            print(f"Warning: Sequence {seqid} not found in genome", file=sys.stderr)
            continue

        seq = genome_seqs[seqid].seq
        py_start = gene_start - 1
        py_end = gene_end

        if strand is Strand.POSITIVE:
            utr5_start = max(0, py_start - upstream_bp[feature.id])
            utr5_end = py_start
            utr3_start = py_end
            utr3_end = min(len(seq), py_end + downstream_bp[feature.id])
        else:
            utr5_start = py_end
            utr5_end = min(len(seq), py_end + upstream_bp[feature.id])
            utr3_start = max(0, py_start - downstream_bp[feature.id])
            utr3_end = py_start

        extracted_regions.extend([
            UTR5Sequence.from_genome(
                seq, seqid, utr5_start, utr5_end, strand, feature.id,
                gene_start=gene_start, gene_end=gene_end
            ),
            CDSSequence.from_genome(
                seq, seqid, py_start, py_end, strand, feature.id,
                gene_start=gene_start, gene_end=gene_end
            ),
            UTR3Sequence.from_genome(
                seq, seqid, utr3_start, utr3_end, strand, feature.id,
                gene_start=gene_start, gene_end=gene_end
            ),
        ])

    return extracted_regions


def get_next_gene_distance(annotation_db, genes, upstream=True):
    """Calculate distance to the next CDS on the SAME STRAND in the specified
    direction for each gene, for use as the full intergenic flanking distance.
    Args: annotation_db: gffutils database object
          genes: list of gffutils.Feature objects for matched genes
          upstream: True for upstream, False for downstream
    Returns: dict {protein_id: distance to next same-strand gene in specified direction}"""

    distances = {}

    for feature in genes:
        if upstream:
            direcion = "<=" if feature.strand == '+' else ">="
        else:
            direcion = ">=" if feature.strand == '+' else "<="

        # AND strand filter ensures we find the nearest CDS on the SAME
        # strand -- the opposite-strand nearest neighbor isn't the right
        # boundary for this gene's own intergenic region.
        query = f"""
                SELECT * FROM features
                WHERE seqid = '{feature.seqid}'
                AND featuretype = 'CDS'
                AND strand = '{feature.strand}'
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
    Args: regions: list of Sequence (or subclass) instances
          output_file: path to output FASTA file"""
    with open(output_file, "w") as f:
        for region in regions:
            f.write(region.to_fasta())


def write_fasta_outputs(regions, utr5_path, cds_path, utr3_path):
    """Write extracted regions to separate FASTA files, split by class.
    Args: regions: list of Sequence (or subclass) instances
          utr5_path: path to output FASTA for 5' UTRs
          cds_path: path to output FASTA for CDS
          utr3_path: path to output FASTA for 3' UTRs"""
    write_fasta_output([r for r in regions if isinstance(r, UTR5Sequence)], utr5_path)
    write_fasta_output([r for r in regions if isinstance(r, CDSSequence)], cds_path)
    write_fasta_output([r for r in regions if isinstance(r, UTR3Sequence)], utr3_path)


def write_bed_output(regions, output_file):
    """Write extracted regions to BED format.
    Args: regions: list of Sequence (or subclass) instances
          output_file: path to output BED file"""
    with open(output_file, 'w') as f:
        f.write("# chrom\tstart\tend\tname\tscore\tstrand\ttype\tgene_start\tgene_end\n")
        for region in regions:
            f.write(region.to_bed_line() + "\n")


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
    Returns: list of NucleotideSequence subclass instances for utr5, cds, and utr3
    """

    gene_coords = get_cds_from_gff(annotation_db, protein_hits)

    print(f"Found coordinates for {len(gene_coords)} CDSs", file=sys.stderr)

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
    Returns: list of CDSSequence instances (seq=None -- no genome FASTA available here;
             use CDSSequence.get_protein_id() + a proteome lookup, or slice the genome
             separately and populate .seq, if the actual sequence is needed)"""

    lower_bound = max(1, motif_start - net_size)
    upper_bound = motif_end + net_size

    cds_features = annotation_db.region(
        seqid=seqid,
        start=lower_bound,
        end=upper_bound,
        featuretype='CDS',
        completely_within=False
    )

    return [
        CDSSequence(
            id=feature.id,
            seq=None,
            seqid=feature.seqid,
            start=feature.start,
            end=feature.end,
            strand=feature.strand,
            attributes=dict(feature.attributes),
        )
        for feature in cds_features
    ]


def extract_cds_protein(protein_id, proteome):
    """Extract protein sequence by ID from proteome dictionary.
    Args: protein_id: protein identifier
          proteome: dictionary of protein sequences (from parse_fasta)
    Returns: protein sequence string or None if not found"""

    if protein_id in proteome:
        return str(proteome[protein_id].seq)

    for seq_id in proteome:
        if protein_id in seq_id or seq_id in protein_id:
            return str(proteome[seq_id].seq)

    return None