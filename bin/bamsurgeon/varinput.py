#!/usr/bin/env python

'''
One input representation for both engines, from either input format.

A VariantRequest is what the legacy whitespace varfiles and a VCF both parse
into. The read-editing engine (snvindel) and the SV engine (addsv) each
consume the fields they care about.

`engine` selects which one handles it, because the type names overlap: an
explicit-allele insertion is an edit to existing reads, while a symbolic
<INS> is a junction to be simulated. Both are called INS.

The legacy readers below are relocated, not rewritten. They are the
compatibility contract for every existing varfile and for
scripts/randomsites.py, so they should be moved and left alone.
'''

import os
import re
import logging

from dataclasses import dataclass, field

import pysam

logger = logging.getLogger(__name__)


SMALL_ENGINE = 'small'
SV_ENGINE = 'sv'

BASES = ('A', 'T', 'C', 'G')

# BIG* was a parallel implementation reached by a size threshold; the
# spellings are kept, the separate implementation is not
BIG_ALIASES = {'BIGDEL': 'DEL', 'BIGINV': 'INV', 'BIGDUP': 'DUP'}


@dataclass
class VariantRequest:
    engine: str          # SMALL_ENGINE or SV_ENGINE
    kind: str            # SNV/INS/DEL for small; DEL/DUP/INV/INS/BND for SV
    chrom: str
    start: int
    end: int
    vaf: float = None

    # small variants
    altbase: str = None
    insseq: str = None

    # structural variants
    ndups: int = 1
    tsdlen: int = 0
    ins_motif: str = None
    insseqfile: str = None
    mate_chrom: str = None
    mate_pos: int = None
    flip_left: bool = False
    flip_right: bool = False

    index: int = 0

    @property
    def is_sv(self):
        return self.engine == SV_ENGINE

    @property
    def is_snv(self):
        return self.engine == SMALL_ENGINE and self.kind == 'SNV'

    def span(self):
        ''' reference interval this request touches, for overlap checks '''
        if self.kind == 'BND':
            return (self.chrom, self.start, self.start + 1)
        return (self.chrom, min(self.start, self.end), max(self.start, self.end) + 1)


# ---------------------------------------------------------------------------
# legacy whitespace varfiles
# ---------------------------------------------------------------------------

def read_snv_varfile(path, maxmuts=0):
    ''' BED-like: chrom start end [vaf] [altbase] '''
    requests = []
    ntried = 0

    with open(path, 'r') as bedfile:
        for line in bedfile:
            if maxmuts and ntried >= maxmuts:
                break
            c = line.strip().split()
            if not c:
                continue

            req = VariantRequest(engine=SMALL_ENGINE, kind='SNV', chrom=c[0],
                                 start=int(c[1]), end=int(c[2]))

            if len(c) > 3:
                req.vaf = float(c[3])
            if len(c) == 5:
                altbase = c[4].upper()
                assert altbase in BASES, "ERROR: ALT " + altbase + " not A, T, C, or G!\n"
                req.altbase = altbase

            requests.append(req)
            ntried += 1

    return requests


def read_indel_varfile(path, maxmuts=0):
    ''' BED-like: chrom start end vaf INS|DEL [seq] '''
    requests = []
    ntried = 0

    with open(path, 'r') as bedfile:
        for line in bedfile:
            if maxmuts and ntried >= maxmuts:
                break
            c = line.strip().split()
            if not c:
                continue

            kind = c[4]
            assert kind in ('INS', 'DEL'), 'indel type must be INS or DEL: %s' % kind

            requests.append(VariantRequest(
                engine=SMALL_ENGINE, kind=kind, chrom=c[0],
                start=int(c[1]), end=int(c[2]), vaf=float(c[3]),
                insseq=c[5] if kind == 'INS' else None))
            ntried += 1

    return requests


def read_sv_varfile(path, default_vaf=1.0, maxmuts=0):
    ''' chrom start end TYPE [type-specific fields] '''
    requests = []
    ntried = 0

    with open(path, 'r') as varfile:
        for line in varfile:
            line = line.strip()
            if line == '' or re.search('^#', line):
                continue
            if maxmuts and ntried >= maxmuts:
                break

            if ';' in line:
                logger.error('compound (";"-chained) mutations are not supported, skipping: %s' % line)
                continue

            requests.append(parse_sv_varfile_line(line, default_vaf))
            ntried += 1

    return requests


def parse_sv_varfile_line(bedline, default_vaf):
    '''
    Split on single whitespace characters rather than runs, because empty
    columns are meaningful: `DUP\\t\\t0.9` means "default copy count, VAF 0.9".
    '''
    fields = re.split(r"[\s]", bedline)

    if len(fields) < 4:
        raise ValueError("Invalid varfile line: %s" % bedline)

    def field(i):
        return fields[i].strip() if len(fields) > i else ''

    def number(i, default):
        v = field(i)
        return type(default)(v) if v != '' else default

    kind = fields[3].upper()

    big = kind in BIG_ALIASES
    if big:
        logger.warning('%s is deprecated, treating as %s' % (kind, BIG_ALIASES[kind]))
        kind = BIG_ALIASES[kind]

    req = VariantRequest(engine=SV_ENGINE, kind=kind, chrom=fields[0],
                         start=int(fields[1]), end=int(fields[2]),
                         vaf=default_vaf)

    if kind == 'INS':
        insspec = field(4)
        assert insspec != '', 'insertion requires a sequence, file, RND or INSLIB: entry'
        if os.path.exists(insspec) or insspec == 'RND' or insspec.startswith('INSLIB:'):
            req.insseqfile = insspec
        else:
            assert re.search('^[ATGCatgc]*$', insspec), "cannot determine SV type: %s" % insspec
            req.insseq = insspec.upper()
        req.tsdlen = number(5, 0)
        motif = field(6)
        if motif:
            assert '^' in motif, 'insertion motif specification requires cut site defined by ^'
            req.ins_motif = motif
        req.vaf = number(7, default_vaf)

    elif kind == 'DUP':
        # BIGDUP never took a copy count; its trailing field is the VAF
        if big:
            req.vaf = number(4, default_vaf)
        else:
            req.ndups = number(4, 1)
            req.vaf = number(5, default_vaf)

    elif kind in ('DEL', 'INV'):
        req.vaf = number(4, default_vaf)

    elif kind in ('TRN', 'BND'):
        req.kind = 'BND'
        req.mate_chrom = field(4)
        req.mate_pos = int(field(5))
        # fields[6] is the mate end; breakends are points, so it is unused
        orient = field(7)
        if orient:
            req.flip_left = orient[0] == '-'
            req.flip_right = orient[1] == '-'
        req.vaf = number(8, default_vaf)

    else:
        raise ValueError("mutation type not one of INS,INV,DEL,DUP,TRN: %s" % kind)

    return req


# ---------------------------------------------------------------------------
# VCF
# ---------------------------------------------------------------------------

def _info(rec, key, default=None):
    '''
    INFO lookup that tolerates undeclared keys. pysam raises
    ValueError('Invalid header') rather than returning None for a key the
    header does not define.
    '''
    try:
        value = rec.info.get(key, default)
    except (ValueError, KeyError):
        return default
    if isinstance(value, (tuple, list)):
        return value[0] if value else default
    return default if value is None else value


def resolve_vaf(rec, default_vaf):
    ''' INFO/VAF, then a single sample's FORMAT/AF, then the CLI default '''
    vaf = _info(rec, 'VAF')
    if vaf is not None:
        return float(vaf)

    if len(rec.samples) == 1:
        for sample in rec.samples.values():
            af = sample.get('AF')
            if af is not None:
                if isinstance(af, (tuple, list)):
                    af = af[0] if af else None
                if af is not None:
                    return float(af)

    return default_vaf


BRACKETS = ('[', ']')


def parse_breakend(alt, chrom, pos):
    '''
    Parse the four bracket forms into the +-/-+ orientation grammar.

        t[p[   base first, '[' -> flip_left False, flip_right False
        t]p]   base first, ']' -> flip_left False, flip_right True
        [p[t   base last,  '[' -> flip_left True,  flip_right False
        ]p]t   base last,  ']' -> flip_left True,  flip_right True

    Returns (mate_chrom, mate_pos, flip_left, flip_right).
    '''
    bracket = None
    for b in BRACKETS:
        if b in alt:
            bracket = b
            break

    if bracket is None:
        raise ValueError('not a breakend ALT: %s' % alt)

    parts = alt.split(bracket)
    if len(parts) != 3:
        raise ValueError('malformed breakend ALT: %s' % alt)

    prefix, locus, suffix = parts

    # the base sits on whichever side is non-empty
    flip_left = (prefix == '' and suffix != '')
    flip_right = (bracket == ']')

    if ':' not in locus:
        raise ValueError('breakend ALT has no mate locus: %s' % alt)

    mate_chrom, mate_pos = locus.rsplit(':', 1)

    return mate_chrom, int(mate_pos), flip_left, flip_right


def read_vcf(path, default_vaf=1.0, maxmuts=0):
    '''
    Read a VCF into requests.

    Explicit REF/ALT goes to the read-editing engine; a symbolic ALT or an
    SVTYPE goes to the SV engine. That is the whole routing rule: the
    representation says which one the author meant.

    Simulation-only directives that have no native VCF home are read from
    BS_-prefixed INFO keys. TSDLEN and NDUPS are read unprefixed because
    that is what makevcf writes, which is what makes output valid input.
    '''
    requests = []
    seen_mates = set()
    ntried = 0

    with pysam.VariantFile(path) as vcf:
        for rec in vcf:
            if maxmuts and ntried >= maxmuts:
                break

            if not rec.alts:
                logger.warning('skipping record with no ALT at %s:%d' % (rec.chrom, rec.pos))
                continue

            alt = rec.alts[0]
            svtype = _info(rec, 'SVTYPE')
            vaf = resolve_vaf(rec, default_vaf)

            req = None

            if svtype == 'BND' or any(b in alt for b in BRACKETS):
                rec_id = rec.id or '%s:%d' % (rec.chrom, rec.pos)
                if rec_id in seen_mates:
                    continue

                mate_chrom, mate_pos, flip_left, flip_right = parse_breakend(
                    alt, rec.chrom, rec.pos)

                # MATEID only dedupes; the mate is synthesised from the ALT
                mateid = _info(rec, 'MATEID')
                if mateid:
                    seen_mates.add(str(mateid))

                req = VariantRequest(
                    engine=SV_ENGINE, kind='BND', chrom=rec.chrom,
                    start=rec.pos, end=rec.pos, vaf=vaf,
                    mate_chrom=mate_chrom, mate_pos=mate_pos,
                    flip_left=flip_left, flip_right=flip_right)

            elif svtype or (alt.startswith('<') and alt.endswith('>')):
                kind = str(svtype or alt.strip('<>')).upper()
                kind = BIG_ALIASES.get(kind, kind)

                if kind not in ('DEL', 'DUP', 'INV', 'INS'):
                    logger.warning('skipping unsupported SVTYPE %s at %s:%d'
                                   % (kind, rec.chrom, rec.pos))
                    continue

                if kind == 'INS':
                    # A varfile gives an interval to place the insertion in and
                    # the midpoint is taken; a VCF gives the exact position it
                    # landed at. Collapsing the interval keeps that midpoint on
                    # POS, so re-running reproduces the same placement.
                    end = rec.pos
                else:
                    end = rec.stop
                    if not end or end <= rec.pos:
                        svlen = _info(rec, 'SVLEN')
                        end = rec.pos + abs(int(svlen)) if svlen else rec.pos + 1

                req = VariantRequest(
                    engine=SV_ENGINE, kind=kind, chrom=rec.chrom,
                    start=rec.pos, end=end, vaf=vaf,
                    ndups=int(_info(rec, 'NDUPS', 1) or 1),
                    tsdlen=int(_info(rec, 'TSDLEN', 0) or 0),
                    ins_motif=_info(rec, 'BS_MOTIF'))

                if kind == 'INS':
                    inslib = _info(rec, 'BS_INSLIB')
                    insseq = _info(rec, 'BS_INSSEQ')
                    if insseq:
                        req.insseq = str(insseq).upper()
                    elif inslib:
                        req.insseqfile = 'INSLIB:' + str(inslib)
                    elif rec.id and rec.id != '.':
                        # older output put the library entry name in ID
                        req.insseqfile = 'INSLIB:' + rec.id

            else:
                ref = rec.ref.upper()
                alt_u = alt.upper()

                if len(ref) == 1 and len(alt_u) == 1:
                    req = VariantRequest(
                        engine=SMALL_ENGINE, kind='SNV', chrom=rec.chrom,
                        start=rec.pos, end=rec.pos, vaf=vaf, altbase=alt_u)

                elif len(alt_u) > len(ref):
                    req = VariantRequest(
                        engine=SMALL_ENGINE, kind='INS', chrom=rec.chrom,
                        start=rec.pos, end=rec.pos + 1, vaf=vaf,
                        insseq=alt_u[len(ref):])

                elif len(ref) > len(alt_u):
                    req = VariantRequest(
                        engine=SMALL_ENGINE, kind='DEL', chrom=rec.chrom,
                        start=rec.pos, end=rec.pos + (len(ref) - len(alt_u)),
                        vaf=vaf)

                else:
                    logger.warning('skipping MNP/complex allele at %s:%d (%s>%s)'
                                   % (rec.chrom, rec.pos, ref, alt))
                    continue

            if req is not None:
                requests.append(req)
                ntried += 1

    return requests


# ---------------------------------------------------------------------------
# input checks
# ---------------------------------------------------------------------------

def sample_read_length(bam_path, fasta_ref=None, limit=1000):
    ''' modal aligned read length over the first `limit` primary reads '''
    counts = {}
    with pysam.AlignmentFile(bam_path, reference_filename=fasta_ref) as bam:
        for i, read in enumerate(bam.fetch(until_eof=True)):
            if i >= limit:
                break
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.query_length:
                counts[read.query_length] = counts.get(read.query_length, 0) + 1

    if not counts:
        return None

    return max(counts.items(), key=lambda kv: (kv[1], kv[0]))[0]


def warn_proximity(clusters, read_length):
    '''
    Warn when sites in different haplotype clusters are close enough to share
    reads.

    Each cluster is spiked independently into its own temp BAM and the results
    are merged, so two clusters touching the same reads overwrite each other
    and one edit is silently lost. --haplosize is what groups them; this
    reports the cases where it is set too low.

    Sites are compared on their requested intervals, not on the position drawn
    from within them, because the draw happens later.
    '''
    flat = []
    for i, cluster in enumerate(clusters):
        for site in cluster:
            flat.append((site.chrom, min(site.start, site.end),
                         max(site.start, site.end), i, site))

    flat.sort(key=lambda t: (t[0], t[1]))

    collisions = 0

    for (chrom_a, _, end_a, ci_a, a), (chrom_b, start_b, _, ci_b, b) in zip(flat, flat[1:]):
        if chrom_a != chrom_b or ci_a == ci_b:
            continue

        gap = start_b - end_a
        if gap >= read_length:
            continue

        collisions += 1
        detail = ''
        if a.vaf is not None and b.vaf is not None and a.vaf != b.vaf:
            detail = ' and their VAFs differ (%s vs %s), so grouping them would collapse both to one' % (a.vaf, b.vaf)

        logger.warning('%s:%d and %s:%d are %dbp apart but in different haplotype '
                       'clusters; they may share reads and overwrite each other%s'
                       % (chrom_a, a.start, chrom_b, b.start, gap, detail))

    if collisions:
        logger.warning('%d site pair(s) closer than one read length (%dbp) are in '
                       'different clusters; consider --haplosize %d or -z auto'
                       % (collisions, read_length, read_length))

    return collisions


def warn_sv_overlap(requests):
    '''
    Warn when a small variant sits inside an SV interval.

    The SV engine excludes every read over its interval and replaces them with
    reads simulated from a reference template, which carries no small variants.
    Anything spiked there first is reverted. Detection only -- composing the
    two is a modelling problem, see doc/design-notes.md.
    '''
    svs = [r for r in requests if r.is_sv]
    small = [r for r in requests if not r.is_sv]

    if not svs or not small:
        return 0

    overlaps = 0

    for s in small:
        s_chrom, s_start, s_end = s.span()
        for v in svs:
            v_chrom, v_start, v_end = v.span()
            if s_chrom == v_chrom and s_start < v_end and v_start < s_end:
                overlaps += 1
                logger.warning('%s at %s:%d falls inside %s %s:%d-%d; SV reads are '
                               'simulated from reference and will revert it'
                               % (s.kind, s.chrom, s.start, v.kind,
                                  v.chrom, v.start, v.end))
                break

    if overlaps:
        logger.warning('%d small variant(s) overlap an SV. Spike SVs first, then '
                       'small variants into the SV output, so the edits survive.'
                       % overlaps)

    return overlaps


# ---------------------------------------------------------------------------
# dispatch
# ---------------------------------------------------------------------------

VCF_SUFFIXES = ('.vcf', '.vcf.gz', '.bcf')


def looks_like_vcf(path):
    if path.lower().endswith(VCF_SUFFIXES):
        return True
    try:
        with open(path, 'rb') as fh:
            return fh.read(16).startswith(b'##fileformat=VCF')
    except OSError:
        return False


def read_variants(path, kind, default_vaf=1.0, maxmuts=0):
    '''
    Read either format. `kind` names the legacy reader to use when the input
    is not a VCF: 'snv', 'indel', 'sv', or 'any' to keep everything, which is
    what the unified entry point wants.
    '''
    if looks_like_vcf(path):
        requests = read_vcf(path, default_vaf=default_vaf, maxmuts=maxmuts)

        if kind == 'any':
            return requests

        # both engines' records are visible here and nowhere else, so this is
        # the only place a single-engine tool can see the cross-engine overlap
        warn_sv_overlap(requests)

        wanted = SV_ENGINE if kind == 'sv' else SMALL_ENGINE
        kept = [r for r in requests if r.engine == wanted]

        if len(kept) != len(requests):
            logger.warning('%s: ignoring %d record(s) not handled by this tool'
                           % (path, len(requests) - len(kept)))

        return kept

    if kind == 'snv':
        return read_snv_varfile(path, maxmuts)
    if kind == 'indel':
        return read_indel_varfile(path, maxmuts)
    if kind == 'sv':
        return read_sv_varfile(path, default_vaf, maxmuts)
    if kind == 'any':
        raise ValueError('a legacy varfile cannot mix variant classes; '
                         'use a VCF for combined input')

    raise ValueError('unknown variant kind: %s' % kind)
