"""
Core sequence classes for DarwinsRNAHunt.

Hierarchy:
    Sequence
     +-- NucleotideSequence
     |    +-- CDSSequence
     |    +-- UTR5Sequence
     |    +-- UTR3Sequence
     |    +-- Motif
     +-- ProteinSequence

Bio.Seq is still used internally to store the actual sequence data; this
module removes the surrounding dict-shaped records and duplicate
header-parsing regexes that used to live in genome_analysis.py.
"""

import re
from dataclasses import dataclass, field
from typing import ClassVar, Optional
from enum import Enum

from Bio.Seq import Seq


class Strand(Enum):
    POSITIVE = "+"
    NEGATIVE = "-"


def _strand(value) -> Strand:
    """Coerce a raw '+'/'-' string or an existing Strand into a Strand."""
    return value if isinstance(value, Strand) else Strand(value)


@dataclass
class Sequence:
    """Base class for any extracted genomic/protein region."""
    id: str
    seq: Seq
    seqid: str
    start: int          # 1-based, inclusive
    end: int             # 1-based, inclusive
    strand: Strand

    region_type: ClassVar[str] = "sequence"

    def __post_init__(self):
        self.strand = _strand(self.strand)
        if not isinstance(self.seq, Seq):
            self.seq = Seq(str(self.seq)) if self.seq is not None else None

    def __len__(self):
        return len(self.seq) if self.seq is not None else 0

    def __str__(self):
        return str(self.seq) if self.seq is not None else ""

    def to_header(self) -> str:
        return f">{self.id}|{self.region_type}|{self.seqid}:{self.start}-{self.end}({self.strand.value})"

    def to_fasta(self) -> str:
        return f"{self.to_header()}\n{str(self)}\n"

    def to_bed_line(self) -> str:
        gene_start = getattr(self, "gene_start", None)
        gene_end = getattr(self, "gene_end", None)
        return (
            f"{self.seqid}\t{self.start - 1}\t{self.end}\t{self.id}|{self.region_type}\t0\t"
            f"{self.strand.value}\t{self.region_type}\t"
            f"{gene_start if gene_start is not None else ''}\t"
            f"{gene_end if gene_end is not None else ''}"
        )


@dataclass
class NucleotideSequence(Sequence):
    gene_start: Optional[int] = None
    gene_end: Optional[int] = None
    attributes: Optional[dict] = None   # raw GFF attributes, when built from a Feature

    region_type: ClassVar[str] = "nucleotide"

    def gap_to(self, other: "NucleotideSequence") -> int:
        """Coordinate gap to another NucleotideSequence. >0 = disjoint gap, <=0 = overlap/contained."""
        return max(self.start - other.end, other.start - self.end)

    def overlaps(self, other: "NucleotideSequence") -> bool:
        return self.gap_to(other) <= 0

    def with_sequence(self, seq) -> "NucleotideSequence":
        """Return a copy of this instance with .seq populated (e.g. after slicing from a genome)."""
        from dataclasses import replace
        return replace(self, seq=seq if isinstance(seq, Seq) else Seq(str(seq)))

    def get_protein_id(self) -> Optional[str]:
        """Pull a protein_id out of GFF attributes (moved from extract_protein_id_from_gff)."""
        if not self.attributes:
            return None
        for key in ['protein_id', 'Name', 'ID', 'Dbxref']:
            if key in self.attributes:
                value = self.attributes[key][0]
                if key == 'Dbxref':
                    for item in value.split(','):
                        if 'protein_id:' in item or 'RefSeq:' in item:
                            return item.split(':')[1]
                else:
                    return value
        return None

    def reverse_complement(self) -> "NucleotideSequence":
        """Return a new instance, revcomp'd, with strand flipped (Bio.Seq doesn't mutate in place)."""
        flipped = Strand.NEGATIVE if self.strand is Strand.POSITIVE else Strand.POSITIVE
        kwargs = dict(id=self.id, seq=self.seq.reverse_complement(), seqid=self.seqid,
                      start=self.start, end=self.end, strand=flipped)
        if hasattr(self, "gene_start"):
            kwargs.update(gene_start=self.gene_start, gene_end=self.gene_end, attributes=self.attributes)
        return type(self)(**kwargs)

    def translate(self, id: Optional[str] = None, assembly_accession: Optional[str] = None) -> "ProteinSequence":
        return ProteinSequence(
            id=id or self.id,
            seq=self.seq.translate(),
            seqid=self.seqid,
            start=self.start,
            end=self.end,
            strand=self.strand,
            assembly_accession=assembly_accession,
        )

    @classmethod
    def from_genome(cls, genome_seq, seqid, py_start, py_end, strand, id, **extra) -> "NucleotideSequence":
        """
        Extract cls's region from a genome sequence.
        py_start/py_end: 0-based python slice coordinates (as gffutils/Biopython use).
        Replaces the free-standing extract_segment function.
        """
        strand = _strand(strand)
        segment = genome_seq[py_start:py_end]
        if strand is Strand.NEGATIVE:
            segment = segment.reverse_complement()
        return cls(id=id, seq=segment, seqid=seqid, start=py_start + 1, end=py_end, strand=strand, **extra)

    @classmethod
    def from_header(cls, header: str, sequence=None) -> Optional["NucleotideSequence"]:
        """
        Parse 'protein_id|type|seqid:start-end(strand)' and dispatch to the
        right subclass based on the type field. Replaces decode_nuc_fasta_record.
        """
        pattern = r'([^|]+)\|([^|]+)\|([^|]+):(\d+)-(\d+)\(([+-])\)'
        match = re.match(pattern, header)
        if not match:
            return None
        region_id, region_type, seqid, start, end, strand = match.groups()
        target_cls = _NUCLEOTIDE_TYPE_MAP.get(region_type, cls)
        return target_cls(id=region_id, seq=sequence, seqid=seqid,
                           start=int(start), end=int(end), strand=strand)


@dataclass
class CDSSequence(NucleotideSequence):
    region_type: ClassVar[str] = "cds"


@dataclass
class UTR5Sequence(NucleotideSequence):
    region_type: ClassVar[str] = "utr5"


@dataclass
class UTR3Sequence(NucleotideSequence):
    region_type: ClassVar[str] = "utr3"


@dataclass
class ProteinSequence(Sequence):
    assembly_accession: Optional[str] = None

    region_type: ClassVar[str] = "Protein"

    def to_header(self) -> str:
        return f">{self.assembly_accession}|{self.id}|{self.region_type}|{self.seqid}:{self.start}-{self.end}({self.strand.value})"

    @classmethod
    def from_cds(cls, nuc: NucleotideSequence, assembly_accession: Optional[str] = None) -> "ProteinSequence":
        return nuc.translate(assembly_accession=assembly_accession)

    @classmethod
    def from_header(cls, header: str, sequence=None) -> Optional["ProteinSequence"]:
        """
        Parse 'assembly|protein_id|type|seqid:start-end(strand)'.
        Replaces decode_prot_fasta_record.
        """
        pattern = r'([^|]+)\|([^|]+)\|([^|]+)\|([^:]+):(\d+)-(\d+)\(([+-])\)'
        match = re.match(pattern, header)
        if not match:
            return None
        assembly, protein_id, _type, seqid, start, end, strand = match.groups()
        return cls(id=protein_id, seq=sequence, seqid=seqid, start=int(start), end=int(end),
                    strand=strand, assembly_accession=assembly)
