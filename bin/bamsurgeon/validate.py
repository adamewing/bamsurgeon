#!/usr/bin/env python

'''
Verify that variants described in a truth VCF are actually present in a BAM.

BAMSurgeon emits a truth VCF alongside every spike-in run, but nothing has ever
checked that the BAM agrees with it. The shell tests end in a `samtools mpileup
| bcftools call -vm` that prints calls and always exits 0.

Everything here is built on pysam directly rather than shelling out to
`samtools mpileup` (which is what mutation.countBaseAtPos does), because the
subprocess route does not work on CRAM without extra plumbing and costs a fork
per position.
'''

import logging
from dataclasses import dataclass, field

import pysam

from bamsurgeon.copynumber import dup_depth_scale

FORMAT = '%(levelname)s %(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


# site statuses
OK = 'OK'                    # found at the requested position and VAF
LOW_VAF = 'LOW_VAF'          # found, but at a materially lower VAF than requested
ABSENT = 'ABSENT'            # covered, but no supporting reads
NO_COVERAGE = 'NO_COVERAGE'  # too little depth to say anything
SHIFTED = 'SHIFTED'          # found, but not where the truth VCF says it is


@dataclass
class Thresholds:
    ''' every cutoff is a knob; defaults are deliberately generous so that a
        pass means something and a failure is worth investigating '''
    min_depth: int = 4           # below this, report NO_COVERAGE rather than ABSENT
    min_alt_reads: int = 2       # below this, report ABSENT
    vaf_ratio_min: float = 0.5   # observed VAF must be >= requested * this
    indel_window: int = None     # default: len(indel) + 1, set per-site
    sv_search_radius: int = 200  # how far to look for a shifted breakend
    sv_exact_tolerance: int = 5  # breakend within this many bp counts as exact
    sv_min_clip: int = 10        # soft-clip length that counts as breakend evidence
    sv_small_bp: int = 50        # at or below this, a DEL/DUP reads as a CIGAR indel
    sv_pass_rate: float = 0.8    # fraction of SV sites that must be OK
    # fraction of a duplication's predicted excess depth that must be present
    dup_depth_min: float = 0.5
    # below this predicted ratio the excess is too small to measure against
    # flank noise, so the check is direction-only
    dup_depth_floor: float = 1.5


@dataclass
class SiteObservation:
    depth: int = 0
    alt_count: int = 0
    alt_fraction: float = 0.0
    status: str = NO_COVERAGE
    # populated for SVs when the breakend is not where it was requested
    observed_pos: int = None
    offset: int = None
    note: str = ''


@dataclass
class SiteReport:
    chrom: str
    pos: int
    kind: str
    requested_vaf: float = None
    observation: SiteObservation = field(default_factory=SiteObservation)

    @property
    def status(self):
        return self.observation.status

    @property
    def passed(self):
        return self.observation.status == OK


def _primary(read):
    return not (read.is_unmapped or read.is_secondary or read.is_supplementary
                or read.is_duplicate or read.is_qcfail)


def _grade(obs, requested_vaf, thresholds):
    ''' shared status logic: depth first, then presence, then VAF '''
    if obs.depth < thresholds.min_depth:
        obs.status = NO_COVERAGE
    elif obs.alt_count < thresholds.min_alt_reads:
        obs.status = ABSENT
    elif requested_vaf is not None and \
            obs.alt_fraction < requested_vaf * thresholds.vaf_ratio_min:
        obs.status = LOW_VAF
    else:
        obs.status = OK
    return obs


def observe_snv(bam, chrom, pos, alt_allele, thresholds, requested_vaf=None):
    ''' pos is 1-based (VCF convention) '''
    obs = SiteObservation()
    start = pos - 1
    alt = alt_allele.upper()

    # truncate=True is essential: without it pysam yields every column spanned
    # by any overlapping read, not just the requested one.
    for col in bam.pileup(chrom, start, start + 1, truncate=True,
                          min_base_quality=0, stepper='nofilter'):
        if col.reference_pos != start:
            continue
        for pread in col.pileups:
            if pread.is_del or pread.is_refskip or pread.query_position is None:
                continue
            if not _primary(pread.alignment):
                continue
            obs.depth += 1
            if pread.alignment.query_sequence[pread.query_position].upper() == alt:
                obs.alt_count += 1

    if obs.depth:
        obs.alt_fraction = obs.alt_count / float(obs.depth)

    return _grade(obs, requested_vaf, thresholds)


def _indel_events(read):
    ''' yield (op, length, refpos) for each I/D in a read's CIGAR.
        refpos is the reference coordinate the event begins at. '''
    if read.cigartuples is None:
        return
    refpos = read.reference_start
    for op, length in read.cigartuples:
        if op in (0, 7, 8):        # M, =, X consume both
            refpos += length
        elif op == 2:              # D
            yield ('D', length, refpos)
            refpos += length
        elif op == 3:              # N
            refpos += length
        elif op == 1:              # I consumes query only
            yield ('I', length, refpos)
        # S(4), H(5), P(6) consume no reference


def observe_indel(bam, chrom, pos, ref_allele, alt_allele, thresholds,
                  requested_vaf=None):
    '''
    Match on (type, length) inside a window rather than on exact position.

    Indels in repeat context are ambiguous under left-alignment: the aligner is
    free to place an equivalent event a few bases either side of where the truth
    VCF puts it. Demanding an exact position match makes this validator report
    failures that are not failures, which is worse than having no validator.

    Insertions get a second detection path. Once the inserted sequence is a
    large fraction of the read length, an aligner stops emitting an I operation
    and soft-clips instead -- a 36bp insertion in 100bp reads comes back as
    60M40S/64S36M, not 36I. So a read also counts as supporting if it simply
    carries the inserted sequence. Matching ~20 exact bases of a known insert is
    specific enough not to fire by chance.
    '''
    obs = SiteObservation()

    inserted = ''
    if len(alt_allele) > len(ref_allele):
        want_op, want_len = 'I', len(alt_allele) - len(ref_allele)
        inserted = alt_allele[len(ref_allele):].upper()
    elif len(ref_allele) > len(alt_allele):
        want_op, want_len = 'D', len(ref_allele) - len(alt_allele)
    else:
        obs.note = 'not an indel'
        return obs

    # long enough to be specific, short enough to survive a sequencing error
    probe = inserted[:20] if len(inserted) >= 8 else ''

    window = thresholds.indel_window
    if window is None:
        window = want_len + 1

    # VCF anchors an indel on the base *before* the event
    expect = pos

    lo = max(0, expect - window - want_len)
    hi = expect + window + want_len + 1

    for read in bam.fetch(chrom, lo, hi):
        if not _primary(read):
            continue
        spans_anchor = read.reference_start <= expect < read.reference_end

        supported = False
        for op, length, refpos in _indel_events(read):
            if op == want_op and length == want_len and abs(refpos - expect) <= window:
                supported = True
                break

        # soft-clipped insertion: the I operation is gone, the sequence is not
        if not supported and probe and read.query_sequence:
            if probe in read.query_sequence.upper():
                supported = True
                obs.note = 'matched by inserted sequence (soft-clipped)'

        if supported:
            obs.alt_count += 1

        # A read that carries the insertion can be clipped so hard that its
        # alignment no longer spans the anchor base. Counting it as support but
        # not as depth yields fractions above 1.0, so the denominator is
        # "covers the site or carries the variant".
        if spans_anchor or supported:
            obs.depth += 1

    if obs.depth:
        obs.alt_fraction = obs.alt_count / float(obs.depth)

    note = obs.note
    obs = _grade(obs, requested_vaf, thresholds)
    if obs.status == OK:
        obs.note = note
    else:
        obs.note = 'looked for %s%d within +/-%dbp of %d%s' % (
            want_op, want_len, window, expect,
            ' (and for the inserted sequence)' if probe else '')
    return obs


def clip_profile(bam, chrom, start, end, min_clip=10):
    '''
    Map reference position -> number of reads whose alignment is soft-clipped
    there. A real breakend shows up as a sharp peak, because every read crossing
    it has to stop aligning at the same base.

    This is the one signal that works identically for DEL, DUP, INV and BND,
    which is what lets observe_sv() treat them uniformly.
    '''
    profile = {}
    start = max(0, start)

    for read in bam.fetch(chrom, start, end):
        if not (_primary(read) or read.is_supplementary):
            continue
        cigar = read.cigartuples
        if not cigar:
            continue
        if cigar[0][0] == 4 and cigar[0][1] >= min_clip:
            profile[read.reference_start] = profile.get(read.reference_start, 0) + 1
        if cigar[-1][0] == 4 and cigar[-1][1] >= min_clip:
            profile[read.reference_end] = profile.get(read.reference_end, 0) + 1

    return profile


def _observe_small_sv(bam, chrom, pos, svlen, svtype, thresholds,
                      requested_vaf=None):
    '''
    Sub-read-length DEL and DUP, detected the way an aligner actually
    represents them: a deletion of N bases is a CIGAR D of length N, and a
    tandem duplication of N bases is a CIGAR I of length N.
    '''
    obs = SiteObservation()
    want_op = 'D' if svtype == 'DEL' else 'I'
    window = svlen + 1

    lo = max(0, pos - window - svlen)
    hi = pos + window + svlen + 1

    for read in bam.fetch(chrom, lo, hi):
        if not _primary(read):
            continue
        if not (read.reference_start <= pos < read.reference_end):
            continue
        obs.depth += 1
        for op, length, refpos in _indel_events(read):
            if op == want_op and length == svlen and abs(refpos - pos) <= window:
                obs.alt_count += 1
                break

    if obs.depth:
        obs.alt_fraction = obs.alt_count / float(obs.depth)

    # Graded on presence, not VAF, for the same reason as the breakend path:
    # only reads that span the whole event can express it, so the supporting
    # fraction sits below the allele fraction by an amount that depends on
    # event length and read length rather than on whether the spike-in worked.
    if obs.depth < thresholds.min_depth:
        obs.status = NO_COVERAGE
    elif obs.alt_count < thresholds.min_alt_reads:
        obs.status = ABSENT
        obs.note = 'short %s: looked for %s%d within +/-%dbp of %d' % (
            svtype, want_op, svlen, window, pos)
    else:
        obs.status = OK
        obs.note = 'short %s matched as CIGAR %s%d' % (svtype, want_op, svlen)

    return obs


def _best_breakend(profile, expect, radius):
    ''' strongest clip peak within radius of expect -> (pos, support) '''
    best_pos, best_support = None, 0
    for pos, support in profile.items():
        if abs(pos - expect) <= radius and support > best_support:
            best_pos, best_support = pos, support
    return best_pos, best_support


def _mean_depth(bam, chrom, start, end):
    start = max(0, start)
    if end <= start:
        return 0.0
    total = sum(sum(col) for col in
                bam.count_coverage(chrom, start, end, quality_threshold=0))
    return total / float(end - start)


def observe_sv(bam, chrom, pos, thresholds, requested_vaf=None,
               orig_bam=None, svtype=None, end=None, tsdlen=0, ndups=None):
    '''
    Locate the breakend by soft-clip peak. If it is not at the requested
    position but is within sv_search_radius, report SHIFTED and the offset --
    that is the status that characterises the assembly-based engine, whose
    coordinates come from an exonerate alignment rather than from the request.

    SVs are deliberately not graded on VAF. The supporting fraction here is
    clipped reads over local depth, which is not an allele fraction: only
    reads that physically span a junction can clip, so even a VAF of 1.0
    produces a fraction well below 1. Comparing it to a requested VAF, as
    this function originally did, made correct spike-ins look like failures.
    The fraction is still reported, and for DEL and DUP the interval depth is
    compared against the original BAM, which is a sound VAF signal and can
    fail a site whose copy number moved the wrong way.
    '''
    obs = SiteObservation()
    expect = pos - 1  # 0-based

    # A short DEL or DUP is an indel-scale event: an aligner emits a CIGAR
    # D or I for it rather than clipping, so the breakend detector below
    # finds nothing to peak on. A 10bp tandem duplication is a 10bp
    # insertion as far as the BAM is concerned.
    svlen = (end - pos) if end else 0
    if svtype in ('DEL', 'DUP') and 0 < svlen <= thresholds.sv_small_bp:
        return _observe_small_sv(bam, chrom, pos, svlen, svtype, thresholds,
                                 requested_vaf)

    radius = thresholds.sv_search_radius
    profile = clip_profile(bam, chrom, expect - radius, expect + radius + 1,
                           min_clip=thresholds.sv_min_clip)

    obs.depth = int(round(_mean_depth(bam, chrom, expect - 50, expect + 50)))

    best_pos, support = _best_breakend(profile, expect, radius)
    obs.alt_count = support
    if obs.depth:
        obs.alt_fraction = support / float(obs.depth)

    if obs.depth < thresholds.min_depth:
        obs.status = NO_COVERAGE
        return obs

    if best_pos is None or support < thresholds.min_alt_reads:
        obs.status = ABSENT
        obs.note = 'no clip peak within %dbp of %d' % (radius, pos)
        return obs

    obs.observed_pos = best_pos + 1  # back to 1-based for reporting
    obs.offset = obs.observed_pos - pos

    # An insertion carrying a target site duplication puts tsdlen bases of
    # reference-matching sequence after POS, so reads stop aligning at the
    # far edge of the TSD rather than at POS itself. Both are the same event.
    tolerated = [0]
    if tsdlen:
        tolerated.append(tsdlen)

    if min(abs(obs.offset - t) for t in tolerated) > thresholds.sv_exact_tolerance:
        obs.status = SHIFTED
        obs.note = 'breakend found at %d (%+d)' % (obs.observed_pos, obs.offset)
        return obs

    obs.status = OK

    # depth corroboration for the two kinds that change copy number; this is
    # the only VAF-sensitive check applied to an SV
    if svtype in ('DEL', 'DUP') and end is not None and orig_bam is not None:
        note, contradicts = _depth_ratio_note(bam, orig_bam, chrom, pos, end,
                                              svtype, thresholds,
                                              requested_vaf, ndups)
        if note:
            obs.note = note
        if contradicts:
            obs.status = LOW_VAF

    return obs


def _depth_ratio_note(bam, orig_bam, chrom, start, end, svtype,
                      thresholds=None, vaf=None, ndups=None):
    '''
    Interval depth in the mutant relative to the original, normalised by the
    same ratio in the flanks. > 1 means sequence was added (DUP), < 1 means it
    was removed (DEL).

    Returns (note, contradicts). A ratio that moved the wrong way is real
    evidence the spike-in did not land, so it can fail a site.

    For a duplication the magnitude is checked too, because VAF and NDUPS
    together predict it: (1 - vaf) + vaf * (ndups + 1). Direction alone passes
    a 4-copy duplication that arrived at 1.02x. The bar is a fraction of the
    predicted *excess* over 1, not of the ratio itself, so it does not get
    easier as the requested duplication gets larger.

    A deletion is left on direction alone. Its prediction, 1 - vaf, is 0 at
    VAF 1.0, so there is no excess to take a fraction of.
    '''
    span = max(1, end - start)
    flank = min(2000, max(200, span // 10))

    try:
        mut_in = _mean_depth(bam, chrom, start, end)
        org_in = _mean_depth(orig_bam, chrom, start, end)
        mut_fl = _mean_depth(bam, chrom, start - flank, start)
        org_fl = _mean_depth(orig_bam, chrom, start - flank, start)
    except ValueError:
        return '', False

    if not (org_in and mut_fl and org_fl):
        return '', False

    normalised = (mut_in / org_in) / (mut_fl / org_fl)
    direction = 'depth ratio %.2f' % normalised

    if svtype == 'DEL' and normalised >= 1.0:
        return direction + ' (expected < 1 for DEL)', True
    if svtype == 'DUP' and normalised <= 1.0:
        return direction + ' (expected > 1 for DUP)', True

    if svtype == 'DUP' and thresholds is not None and vaf is not None:
        predicted = dup_depth_scale(vaf, ndups)

        if predicted >= thresholds.dup_depth_floor:
            floor = 1.0 + (predicted - 1.0) * thresholds.dup_depth_min
            direction += ' of %.2f predicted' % predicted
            if normalised < floor:
                return direction + ' (below %.2f)' % floor, True

    return direction, False


def classify(rec):
    '''
    Decide what kind of variant a VCF record describes.

    BAMSurgeon's own writers are the primary input: makevcf.write_vcf_snv emits
    bare REF/ALT, write_vcf_indel emits bare REF/ALT with differing lengths, and
    write_vcf_sv emits SVTYPE plus a symbolic or breakend ALT.
    '''
    svtype = info_get(rec, 'SVTYPE')
    if svtype:
        return str(svtype)

    alt = rec.alts[0] if rec.alts else ''

    if alt.startswith('<') and alt.endswith('>'):
        return alt.strip('<>')

    if '[' in alt or ']' in alt:
        return 'BND'

    if len(rec.ref) == 1 and len(alt) == 1:
        return 'SNV'

    return 'INDEL'


def info_get(rec, key, default=None):
    '''
    INFO lookup that tolerates undeclared keys.

    pysam raises ValueError('Invalid header') for a key the VCF header does
    not define, rather than returning None, so a truth VCF written by an
    older bamsurgeon blows up on fields added later.
    '''
    try:
        value = rec.info.get(key, default)
    except (ValueError, KeyError):
        return default

    if isinstance(value, (tuple, list)):
        return value[0] if value else default

    return default if value is None else value


def requested_vaf(rec):
    ''' bamsurgeon writes VAF= into INFO on every record it emits '''
    vaf = info_get(rec, 'VAF')
    if vaf is None:
        return None
    if isinstance(vaf, (tuple, list)):
        vaf = vaf[0]
    try:
        return float(vaf)
    except (TypeError, ValueError):
        return None


def validate_vcf(truth_vcf, mutant_bam, ref_fasta, orig_bam=None,
                 thresholds=None):
    ''' returns a list of SiteReport, one per record in truth_vcf '''
    thresholds = thresholds or Thresholds()

    reports = []

    vcf = pysam.VariantFile(truth_vcf)
    bam = pysam.AlignmentFile(mutant_bam, reference_filename=ref_fasta)
    orig = pysam.AlignmentFile(orig_bam, reference_filename=ref_fasta) \
        if orig_bam else None

    try:
        for rec in vcf:
            kind = classify(rec)
            vaf = requested_vaf(rec)
            report = SiteReport(chrom=rec.chrom, pos=rec.pos, kind=kind,
                                requested_vaf=vaf)

            alt = rec.alts[0] if rec.alts else ''

            try:
                if kind == 'SNV':
                    report.observation = observe_snv(
                        bam, rec.chrom, rec.pos, alt, thresholds, vaf)
                elif kind == 'INDEL':
                    report.observation = observe_indel(
                        bam, rec.chrom, rec.pos, rec.ref, alt, thresholds, vaf)
                else:
                    # htslib folds INFO/END into the record's stop, so
                    # info['END'] reads back as None even when the text
                    # carries it. rec.stop is the reliable accessor.
                    end = rec.stop
                    if not end or end <= rec.pos:
                        end = info_get(rec, 'END')

                    tsdlen = info_get(rec, 'TSDLEN', 0) or 0

                    report.observation = observe_sv(
                        bam, rec.chrom, rec.pos, thresholds, vaf,
                        orig_bam=orig, svtype=kind, end=end, tsdlen=int(tsdlen),
                        ndups=info_get(rec, 'NDUPS'))
            except ValueError as e:
                # contig absent from the BAM, or a coordinate off the end of it
                report.observation = SiteObservation(
                    status=NO_COVERAGE, note=str(e))

            reports.append(report)
    finally:
        vcf.close()
        bam.close()
        if orig is not None:
            orig.close()

    return reports


def has_bs_tag_support(mutant_bam, ref_fasta, chrom, pos, alt_allele):
    '''
    For runs made with --tagreads, every altered read carries a BS tag. Same
    check scripts/remove_non_BS.py performs, reused here so that a --tagreads
    run can assert the tagging actually happened.
    '''
    tagged = untagged = 0
    bam = pysam.AlignmentFile(mutant_bam, reference_filename=ref_fasta)
    try:
        start = pos - 1
        for col in bam.pileup(chrom, start, start + 1, truncate=True,
                              min_base_quality=0, stepper='nofilter'):
            if col.reference_pos != start:
                continue
            for pread in col.pileups:
                if pread.is_del or pread.is_refskip or pread.query_position is None:
                    continue
                read = pread.alignment
                if not _primary(read):
                    continue
                base = read.query_sequence[pread.query_position].upper()
                if base != alt_allele.upper():
                    continue
                if read.has_tag('BS'):
                    tagged += 1
                else:
                    untagged += 1
    finally:
        bam.close()

    return tagged, untagged


def summarise(reports):
    ''' status -> count, plus overall pass rate '''
    counts = {}
    for r in reports:
        counts[r.status] = counts.get(r.status, 0) + 1
    total = len(reports)
    passed = counts.get(OK, 0)
    return counts, passed, total
