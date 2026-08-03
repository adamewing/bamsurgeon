#!/usr/bin/env python

'''
Assert that a BAMSurgeon spike-in actually landed.

Reads the truth VCF that addsnv/addindel/addsv emit, looks for each variant in
the output BAM, and exits nonzero if too few were found. This is the assertion
the shell tests never had -- they end in `samtools mpileup | bcftools call -vm`,
which prints calls and always exits 0.

    validate_spikein.py --vcf truth.vcf --bam mutant.bam \\
        --orig-bam input.bam --ref ref.fa [--tsv report.tsv]
'''

import argparse
import os
import sys

# allow running straight out of a checkout, without installing
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'bin'))

from bamsurgeon.validate import (  # noqa: E402
    OK, LOW_VAF, ABSENT, NO_COVERAGE, SHIFTED,
    Thresholds, validate_vcf, summarise, has_bs_tag_support,
)

SV_KINDS = ('DEL', 'DUP', 'INV', 'INS', 'BND', 'TRN')

HEADER = ('chrom', 'pos', 'kind', 'req_vaf', 'depth', 'alt', 'obs_vaf',
          'status', 'note')


def rows(reports):
    for r in reports:
        o = r.observation
        yield (
            r.chrom,
            str(r.pos),
            r.kind,
            '%.3f' % r.requested_vaf if r.requested_vaf is not None else '.',
            str(o.depth),
            str(o.alt_count),
            '%.3f' % o.alt_fraction,
            o.status,
            o.note or '',
        )


def print_table(reports, stream=sys.stdout):
    table = [HEADER] + list(rows(reports))
    widths = [max(len(row[i]) for row in table) for i in range(len(HEADER))]
    for row in table:
        stream.write('  '.join(cell.ljust(widths[i])
                               for i, cell in enumerate(row)).rstrip() + '\n')


def write_tsv(reports, path):
    with open(path, 'w') as fh:
        fh.write('\t'.join(HEADER) + '\n')
        for row in rows(reports):
            fh.write('\t'.join(row) + '\n')


def main():
    p = argparse.ArgumentParser(
        description='verify a BAMSurgeon spike-in against its truth VCF')
    p.add_argument('--vcf', required=True, help='truth VCF emitted by bamsurgeon')
    p.add_argument('--bam', required=True, help='BAM/CRAM produced by the spike-in')
    p.add_argument('--ref', required=True, help='reference FASTA (samtools faidx indexed)')
    p.add_argument('--orig-bam', default=None,
                   help='the input BAM, used to baseline SV interval depth')
    p.add_argument('--tsv', default=None, help='also write the per-site table here')

    p.add_argument('--min-depth', type=int, default=4,
                   help='below this depth a site reports NO_COVERAGE (default 4)')
    p.add_argument('--min-alt-reads', type=int, default=2,
                   help='below this many supporting reads a site is ABSENT (default 2)')
    p.add_argument('--vaf-ratio-min', type=float, default=0.5,
                   help='observed VAF must be at least requested*this (default 0.5)')
    p.add_argument('--indel-window', type=int, default=None,
                   help='bp either side of POS to accept an indel (default: len+1)')
    p.add_argument('--sv-search-radius', type=int, default=200,
                   help='how far to look for a shifted breakend (default 200)')
    p.add_argument('--sv-exact-tolerance', type=int, default=5,
                   help='breakend within this many bp counts as exact (default 5)')
    p.add_argument('--sv-min-clip', type=int, default=10,
                   help='soft-clip length counting as breakend evidence (default 10)')

    p.add_argument('--pass-rate', type=float, default=1.0,
                   help='fraction of SNV/indel sites that must be OK (default 1.0)')
    p.add_argument('--sv-pass-rate', type=float, default=0.8,
                   help='fraction of SV sites that must be OK (default 0.8)')
    p.add_argument('--require-bs-tag', action='store_true', default=False,
                   help='for --tagreads runs, assert alt reads carry the BS tag')

    args = p.parse_args()

    for path, what in ((args.vcf, 'vcf'), (args.bam, 'bam'), (args.ref, 'ref')):
        if not os.path.exists(path):
            sys.exit('no such %s: %s' % (what, path))

    thresholds = Thresholds(
        min_depth=args.min_depth,
        min_alt_reads=args.min_alt_reads,
        vaf_ratio_min=args.vaf_ratio_min,
        indel_window=args.indel_window,
        sv_search_radius=args.sv_search_radius,
        sv_exact_tolerance=args.sv_exact_tolerance,
        sv_min_clip=args.sv_min_clip,
        sv_pass_rate=args.sv_pass_rate,
    )

    reports = validate_vcf(args.vcf, args.bam, args.ref,
                           orig_bam=args.orig_bam, thresholds=thresholds)

    if not reports:
        sys.exit('no records in %s -- nothing to validate' % args.vcf)

    print_table(reports)
    if args.tsv:
        write_tsv(reports, args.tsv)

    counts, passed, total = summarise(reports)

    print('')
    for status in (OK, LOW_VAF, SHIFTED, ABSENT, NO_COVERAGE):
        if counts.get(status):
            print('%-12s %d' % (status, counts[status]))
    print('%-12s %d/%d (%.1f%%)' % ('pass rate', passed, total,
                                    100.0 * passed / total))

    # SVs and small variants are held to different bars: SV breakend detection
    # is inherently noisier, and the SV tests are aggregate by design.
    sv = [r for r in reports if r.kind in SV_KINDS]
    small = [r for r in reports if r.kind not in SV_KINDS]

    failures = []

    if small:
        rate = sum(1 for r in small if r.passed) / float(len(small))
        print('%-12s %.3f (%d sites, need %.3f)'
              % ('snv/indel', rate, len(small), args.pass_rate))
        if rate < args.pass_rate:
            failures.append('snv/indel pass rate %.3f < %.3f'
                            % (rate, args.pass_rate))

    if sv:
        rate = sum(1 for r in sv if r.passed) / float(len(sv))
        print('%-12s %.3f (%d sites, need %.3f)'
              % ('sv', rate, len(sv), args.sv_pass_rate))
        if rate < args.sv_pass_rate:
            failures.append('sv pass rate %.3f < %.3f' % (rate, args.sv_pass_rate))

    if args.require_bs_tag:
        tagged = untagged = 0
        for r in reports:
            if r.kind == 'SNV' and r.passed:
                t, u = has_bs_tag_support(args.bam, args.ref, r.chrom, r.pos,
                                          _alt_of(args.vcf, r))
                tagged += t
                untagged += u
        print('%-12s %d tagged, %d untagged' % ('BS tag', tagged, untagged))
        if untagged:
            failures.append('%d alt-supporting reads lack the BS tag' % untagged)

    if failures:
        print('')
        for f in failures:
            print('FAIL: %s' % f)
        return 1

    return 0


def _alt_of(vcf_path, report):
    ''' re-read the ALT for one site; only used on the --require-bs-tag path '''
    import pysam
    with pysam.VariantFile(vcf_path) as vcf:
        for rec in vcf.fetch(report.chrom, report.pos - 1, report.pos):
            if rec.pos == report.pos and rec.alts:
                return rec.alts[0]
    return 'N'


if __name__ == '__main__':
    sys.exit(main())
