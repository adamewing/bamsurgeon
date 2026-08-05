'''
Local copy number, for scaling the depth thresholds.

--maxdepth is an absolute cap, and duplicated sequence legitimately carries
more reads than the rest of the genome. That matters most in a unified run,
which applies SVs first: a small variant inside a 3-copy duplication is
evaluated against depth that is ~3.4x higher than the input BAM's, and at a
600x starting depth it crosses the 2000x default and is silently dropped with
a threshold message. The same applies to a region a --cnvfile already
describes as amplified.

Only --maxdepth is scaled. Scaling --mindepth would make it stricter over
duplicated sequence and could newly drop sites that pass today.
'''

import logging

import pysam

logger = logging.getLogger(__name__)

DIPLOID = 2.0


def dup_depth_scale(vaf, ndups):
    '''
    Expected depth over a tandem duplication, relative to its flanks.

    A fraction `vaf` of the reads carry `ndups + 1` copies and the rest carry
    one, so the mean is (1 - vaf) + vaf * (ndups + 1). Measured against the
    engine's own output: ndups=3 at VAF 0.9 predicts 3.7 and reads 3.42,
    ndups=1 at VAF 1.0 predicts 2.0 and reads 1.95.
    '''
    vaf = 1.0 if vaf is None else float(vaf)
    ndups = 1 if not ndups else int(ndups)

    return (1.0 - vaf) + vaf * (ndups + 1)


def depth_scale(chrom, start, end, cnvfile=None, intervals=()):
    '''
    Largest depth multiplier that applies anywhere over [start, end).

    `intervals` are (chrom, start, end, scale) tuples, which is how a unified
    run passes down the duplications it applied earlier in the same run.
    `cnvfile` is the tabix-indexed absolute copy number list --cnvfile already
    takes, where 2 means no change.

    The maximum rather than the product: overlapping records describe the same
    locus, so multiplying them would compound one amplification into several.
    '''
    scale = 1.0

    for ivl_chrom, ivl_start, ivl_end, ivl_scale in intervals:
        if ivl_chrom == chrom and ivl_start < end and start < ivl_end:
            scale = max(scale, float(ivl_scale))

    if cnvfile:
        try:
            cnv = pysam.Tabixfile(cnvfile, 'r')
        except (IOError, OSError, ValueError):
            logger.warning('could not read copy number from %s' % cnvfile)
            return scale

        if chrom in cnv.contigs:
            for cnregion in cnv.fetch(chrom, max(0, start), end):
                cn = float(cnregion.strip().split()[3])  # chrom, start, end, CN
                if cn > 0.0:
                    scale = max(scale, cn / DIPLOID)

    return scale


def dup_intervals(records):
    '''
    (chrom, start, end, scale) for every duplication in `records`.

    MutationRecord carries 1-based POS and inclusive END; these come back
    0-based half-open, to match the coordinates a spike-in site is given in.
    '''
    return [(r.chrom, r.pos - 1, r.end, dup_depth_scale(r.vaf, r.ndups))
            for r in records if r.kind == 'DUP']
