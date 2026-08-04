#!/usr/bin/env python

'''
What a completed SV spike-in reports back to the parent process.

addsv runs each mutation in a ProcessPoolExecutor, so anything returned has
to pickle. This is a frozen dataclass of primitives for that reason -- no
pysam handles, no open files.

Previously the SV path communicated results by writing tab-delimited lines
into per-mutation log files, which makevcf then re-parsed with positional
indexing that differed per mutation type. makemut already returned a mutinfo
dict that main() unpacked and discarded; these records take that return path
and make it the real one.
'''

from dataclasses import dataclass, field


# symbolic SV interval types carry START/END; BND carries a mate
INTERVAL_KINDS = ('DEL', 'DUP', 'INV', 'INS')

# small variants are written with explicit REF/ALT rather than a symbolic ALT
SMALL_KINDS = ('SNV', 'INDEL')


@dataclass(frozen=True)
class MutationRecord:
    kind: str          # DEL, DUP, INV, INS, BND
    chrom: str
    pos: int           # 1-based VCF POS; equals the requested start
    end: int = 0       # 1-based inclusive END for interval kinds
    svlen: int = 0
    vaf: float = 1.0

    # BND only
    mate_chrom: str = None
    mate_pos: int = None
    flip_left: bool = False
    flip_right: bool = False

    # INS only
    ins_id: str = None      # insertion library entry, when one was used
    insseq: str = None      # literal inserted sequence, when given inline
    ins_motif: str = None   # preferred cut site, NNN^NNN
    tsdlen: int = 0

    # DUP only
    ndups: int = 1

    # small variants (SNV/INDEL) carry explicit alleles instead of a symbolic
    # ALT, and report the coverage the spike-in achieved
    ref_allele: str = None
    alt_allele: str = None
    dpr: float = None

    # provenance
    mutid: str = ''
    donor_interior: bool = False   # DUP interior came from --donorbam

    @property
    def is_small(self):
        return self.ref_allele is not None and self.alt_allele is not None

    def key(self):
        ''' stable identity, used to derive breakend IDs '''
        return '|'.join(str(x) for x in (
            self.kind, self.chrom, self.pos, self.end,
            self.mate_chrom, self.mate_pos))


@dataclass
class MutationResult:
    ''' what makemut hands back: a temp BAM, reads to drop, and the records '''
    bamfile: str = None
    excludefile: str = None
    records: list = field(default_factory=list)

    def ok(self):
        return self.bamfile is not None and self.excludefile is not None
