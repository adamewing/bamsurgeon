#!/usr/bin/env python

'''
Build the mutated haplotype for an SV directly from reference slices.

The assembly-based engine had to discover where its contig sat before it could
place a breakpoint, which is what align()/trim_contig()/locate_contig_pos()
existed for. Here the interval is chosen rather than discovered, so every
offset is exact by construction and the reported coordinates are the requested
ones.

Coordinates follow the convention the varfile and the emitted VCF already use:
start and end are 1-based positions, VCF POS equals the requested start, and
the affected bases are start+1..end inclusive (start is the anchor base, as
VCF requires for symbolic ALTs).

The sequence edits themselves are bamsurgeon.mutableseq operations, unchanged.
They always worked on plain strings; only the offsets were ever in doubt.
'''

import logging
import re

import bamsurgeon.mutableseq as ms

logger = logging.getLogger(__name__)


class TemplateError(Exception):
    ''' raised when a template cannot be built (N content, off the contig) '''


class Template:
    '''
    seq      the mutated haplotype to hand to wgsim
    exclude  reference intervals whose original reads this replaces
    '''
    def __init__(self, seq, exclude, event_pos=None, event_len=0):
        self.seq = seq
        self.exclude = exclude
        # where the event actually landed, for kinds whose reported position
        # is not simply the requested start (an insertion is placed within
        # the requested interval, not at its edge)
        self.event_pos = event_pos
        self.event_len = event_len

    @property
    def exclude_length(self):
        return sum(e - s for _, s, e in self.exclude)

    def reads_ratio(self):
        '''
        How many reads to simulate, relative to how many were removed.

        Preserving read density across the swap means
            nsim / len(template) == nexcluded / len(replaced reference)
        A deletion shrinks this below 1, a duplication pushes it above 1. The
        old code applied a ratio on the non-translocation path only and left
        the translocation path unscaled.
        '''
        if self.exclude_length <= 0:
            return 0.0
        return len(self.seq) / float(self.exclude_length)


def _fetch(reffile, chrom, start_1, end_1):
    ''' 1-based inclusive fetch, clamped to the contig '''
    length = reffile.get_reference_length(chrom)
    start_1 = max(1, start_1)
    end_1 = min(length, end_1)
    if end_1 <= start_1:
        raise TemplateError('empty interval %s:%d-%d' % (chrom, start_1, end_1))
    return reffile.fetch(chrom, start_1 - 1, end_1).upper(), start_1


def _check_n(seq, mutid, allow_n):
    '''
    Reference N content is known before any work happens, so this is a
    pre-check rather than something discovered after assembling.

    This is also why it now behaves predictably: velvet inserts N runs when it
    scaffolds a paired library across a coverage gap, so the old check could
    fire on intervals whose reference contains no N at all.
    '''
    if 'N' not in seq:
        return seq
    if allow_n:
        logger.warning("%s reference slice has ambiguous bases (N), replaced with 'A'" % mutid)
        return re.sub('N', 'A', seq)
    raise TemplateError('reference slice contains N (use --allowN to override)')


def _offset(pos_1, slice_start_1):
    ''' 0-based index within the slice of the base at 1-based pos_1 '''
    return pos_1 - slice_start_1


def build_interval(kind, reffile, chrom, start, end, pad, mutid='null',
                   allow_n=False, insseq='', tsdlen=0, ndups=1, ins_motif=None,
                   maxlibsize=600):
    '''
    DEL, INV, DUP and INS all reduce to one reference slice plus one edit.

    The slice spans start-pad .. end+pad so that simulated pairs crossing each
    junction are fully contained in the template.
    '''
    seq, slice_start = _fetch(reffile, chrom, start - pad, end + pad)
    seq = _check_n(seq, mutid, allow_n)

    mutseq = ms.MutableSeq(seq)

    event_pos = start
    event_len = end - start

    a = _offset(start, slice_start) + 1   # first affected base
    b = _offset(end, slice_start) + 1     # one past the last affected base

    # an insertion is a point event, so start == end is legitimate there;
    # the other kinds need a non-empty interval to act on
    min_span = 0 if kind == 'INS' else 1

    if a < 1 or b > len(seq) or b - a < min_span:
        raise TemplateError('%s: interval does not fit its slice (a=%d b=%d len=%d)'
                            % (mutid, a, b, len(seq)))

    if kind == 'DEL':
        mutseq.deletion(a, b)

    elif kind == 'INV':
        mutseq.inversion(a, b)

    elif kind == 'DUP':
        mutseq.duplication(a, b, ndups)

    elif kind == 'INS':
        # the varfile gives an interval for insertions; take its midpoint, as
        # the assembly-based path did once it had mapped the interval back
        loc = _offset((start + end) // 2, slice_start) + 1
        if ins_motif is not None:
            # find_site works on any string; it never needed a contig
            trim = int(maxlibsize)
            found = mutseq.find_site(ins_motif, left_trim=trim, right_trim=trim)
            if trim <= found <= mutseq.length() - trim:
                loc = found
            else:
                logger.info('%s no cutsite found, using requested position' % mutid)
        mutseq.insertion(loc, insseq, tsdlen)

        # report where it went, not the interval it was allowed to go in
        event_pos = slice_start + loc - 1
        event_len = len(insseq) + tsdlen

    else:
        raise TemplateError('not an interval mutation: %s' % kind)

    return Template(mutseq.seq, [(chrom, start - pad, end + pad)],
                    event_pos=event_pos, event_len=event_len)


def build_dup_junction(reffile, chrom, start, end, pad, mutid='null',
                       allow_n=False):
    '''
    Tandem duplication as a junction only: the end of one copy joined to the
    start of the next. Used when --donorbam supplies the interior instead of
    simulating it, which keeps the template small for large duplications and
    gives the interior a real error profile.
    '''
    left, _ = _fetch(reffile, chrom, end - pad + 1, end)
    right, _ = _fetch(reffile, chrom, start, start + pad - 1)

    seq = _check_n(left + right, mutid, allow_n)

    return Template(seq, [(chrom, end - pad + 1, end),
                          (chrom, start, start + pad - 1)])


def build_breakend(reffile, chrom1, pos1, chrom2, pos2, pad, mutid='null',
                   allow_n=False, flip_left=False, flip_right=False):
    '''
    Fusion of two reference slices. Only the retained half of each side is
    excluded: the discarded halves still exist on the reference allele and
    their reads must survive.
    '''
    left_seq, left_start = _fetch(reffile, chrom1, pos1 - pad, pos1 + pad)
    right_seq, right_start = _fetch(reffile, chrom2, pos2 - pad, pos2 + pad)

    left_seq = _check_n(left_seq, mutid, allow_n)
    right_seq = _check_n(right_seq, mutid, allow_n)

    left = ms.MutableSeq(left_seq)
    right = ms.MutableSeq(right_seq)

    loc1 = _offset(pos1, left_start) + 1   # keep up to and including pos1
    loc2 = _offset(pos2, right_start)      # resume at pos2

    left.fusion(loc1, right, loc2, flip1=flip_left, flip2=flip_right)

    if flip_left:
        excl_left = (chrom1, pos1, pos1 + pad)
    else:
        excl_left = (chrom1, pos1 - pad, pos1)

    if flip_right:
        excl_right = (chrom2, pos2 - pad, pos2)
    else:
        excl_right = (chrom2, pos2, pos2 + pad)

    return Template(left.seq, [excl_left, excl_right])


def sample_read_length(bamfile, chrom, start, end, limit=1000):
    '''
    Modal aligned read length over the region.

    Replaces the old maxqlen, which was the longest BAM quality string in the
    region, carried through the contig object purely to be len()-ed. The mode
    is a better read length for long-read data, where the maximum is an
    outlier rather than a representative.
    '''
    counts = {}
    for i, read in enumerate(bamfile.fetch(chrom, start, end)):
        if i >= limit:
            break
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.query_length:
            counts[read.query_length] = counts.get(read.query_length, 0) + 1

    if not counts:
        return None

    return max(counts.items(), key=lambda kv: (kv[1], kv[0]))[0]
