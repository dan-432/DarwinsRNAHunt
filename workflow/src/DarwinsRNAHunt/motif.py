from enum import Enum
from Bio.Seq import Seq
from DarwinsRNAHunt.sequence import NucleotideSequence, Strand
from dataclasses import dataclass
from typing import ClassVar


POSITION_COLUMNS = [
    "protein_id",
    "position",                 # upstream / downstream / contained
    "rank",                     # 0 for contained; 1..N within subset otherwise
    "bp_distance_from_motif",
    "cds_strand",
    "motif_strand",
    "same_strand",
    "seqid",
    "cds_start",
    "cds_end",
    "motif_fully_intergenic",
]

def _parse_bool(value):
    """Parse pandas-style 'True'/'False' strings (and tolerate blanks/odd casing)."""
    if value is None:
        return None
    value = value.strip()
    if value == "":
        return None
    return value.strip().lower() == "true"


def _parse_float(value):
    if value is None:
        return None
    value = value.strip()
    if value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None

def _parse_header(value):
    """removed > from the start if there is one"""
    value = value.strip()

    value = value.lstrip(">")

    return value


class Position(Enum):
    UPSTREAM = "upstream"
    DOWNSTREAM = "downstream"
    CONTAINED = "contained"

class CDSPosition:
    """Minimal: describes one CDS's position relative to *a* motif, without
    holding a reference to which one."""
    def __init__(self, record, distance, position, same_strand, rank):
        self.cds = record
        self.distance = distance
        self.position = position
        self.same_strand = same_strand
        self.rank = rank
        self.protein = None


@dataclass
class Motif(NucleotideSequence):
    fully_intergenic: bool = False

    region_type: ClassVar[str] = "motif"

    @classmethod
    def from_hit(cls, id: str, raw_start: int, raw_end: int, sequence: str, seq_id: str) -> "Motif":
        """
        cmsearch convention: raw_start > raw_end signals a minus-strand hit.
        Normalizes coordinates and revcomps the sequence accordingly.
        """
        strand = Strand.NEGATIVE if raw_start > raw_end else Strand.POSITIVE
        start, end = min(raw_start, raw_end), max(raw_start, raw_end)
        seq = Seq(sequence)
        if strand is Strand.NEGATIVE:
            seq = seq.reverse_complement()
        return cls(id=id, seq=seq, seqid=seq_id, start=start, end=end, strand=strand)

    def get_cds_positions_relative_to(self, cds_records):
        """
        ========================= PLUS (+) STRAND =========================
        
        5'                                                              3'
        ----[=====CDS=====]----((motif))----[=====CDS=====]---- - - ->
            s             e     s      e     s            e
        
        if CDS.strand == '+':
            if CDS.e < motif.s:   -> downstream
            if CDS.s > motif.e:   -> upstream
            else:                 -> contained
        
        
        ========================= MINUS (-) STRAND =========================
        
        3'                                                              5'
        <- - -   ----[=====CDS=====]---((motif))---[=====CDS=====]----
                     s             e    s     e    s             e
        
        if CDS.strand == '-':
            if CDS.e < motif.s:   -> upstream
            if CDS.s > motif.e:   -> downstream
            else:                 -> contained
        
        
        Args:
            cds_records: list of NucleotideSequence records of surrounding CDSs
        Returns:
            list of CDSPosition, one per input cds_record, each with
            .position/.same_strand/.rank set (rank=0 for contained)
        """

        positions = []
        self.fully_intergenic = True

        for record in cds_records:
            same_strand = (record.strand == self.strand)
            distance = self.gap_to(record)
            position = CDSPosition(record, distance, Position.CONTAINED, same_strand, 0)

            if record.strand is Strand.POSITIVE:
                if record.end < self.start:
                    position.position = Position.DOWNSTREAM
                elif record.start > self.end:
                    position.position = Position.UPSTREAM
                else:
                    position.position = Position.CONTAINED
                    self.fully_intergenic = False
            else:
                if record.end < self.start:
                    position.position = Position.UPSTREAM
                elif record.start > self.end:
                    position.position = Position.DOWNSTREAM
                else:
                    position.position = Position.CONTAINED
                    self.fully_intergenic = False

            positions.append(position)

        for same_strand_flag, pos_type in [(True, Position.UPSTREAM), (True, Position.DOWNSTREAM),
                                            (False, Position.UPSTREAM), (False, Position.DOWNSTREAM)]:
            subset = [p for p in positions if p.same_strand == same_strand_flag and p.position == pos_type]
            subset.sort(key=lambda p: p.distance)
            for rank, p in enumerate(subset, start=1):
                p.rank = rank

        self.cds_positions = positions

        return positions

    def to_rows(self) -> list[dict]:
        """One row per CDSPosition with a resolved protein. Call
        get_cds_positions_relative_to() first."""
        return [
            PositionRecord(
                protein_id=p.protein.to_header(),
                position=p.position,
                rank=p.rank,
                distance=p.distance,
                cds_strand=p.cds.strand,
                motif_strand=self.strand,
                same_strand=p.same_strand,
                seqid=p.cds.seqid,
                cds_start=p.cds.start,
                cds_end=p.cds.end,
                motif_fully_intergenic=self.fully_intergenic
            ).to_row()
            for p in self.cds_positions if p.protein is not None
        ]


@dataclass
class PositionRecord:
    """A single row from the position TSV, read back for analysis.
    Deliberately NOT a CDSPosition or Motif -- those model the live
    extraction pipeline; this models a stored, flattened record."""
    protein_id: str
    position: Position
    rank: int
    distance: int
    cds_strand: Strand
    motif_strand: Strand
    same_strand: bool
    seqid: str
    cds_start: int
    cds_end: int
    motif_fully_intergenic: bool

    @classmethod
    def from_row(cls, row: dict) -> "PositionRecord":
        return cls(
            protein_id=_parse_header(row["protein_id"]),
            position=Position(row["position"]),
            rank=int(row["rank"]),
            distance=int(row["bp_distance_from_motif"]),
            cds_strand=Strand(row["cds_strand"]),
            motif_strand=Strand(row["motif_strand"]),
            same_strand=_parse_bool(row["same_strand"]),
            seqid=row["seqid"],
            cds_start=int(row["cds_start"]),
            cds_end=int(row["cds_end"]),
            motif_fully_intergenic=_parse_bool(row["motif_fully_intergenic"])
        )

    def to_row(self):
        return {
            "protein_id": self.protein_id,
            "position": self.position.value,
            "rank": self.rank,
            "bp_distance_from_motif": self.distance,
            "cds_strand": self.cds_strand.value,
            "motif_strand": self.motif_strand.value,
            "same_strand": self.same_strand,
            "seqid": self.seqid,
            "cds_start": self.cds_start,
            "cds_end": self.cds_end,
            "motif_fully_intergenic": self.motif_fully_intergenic
        }