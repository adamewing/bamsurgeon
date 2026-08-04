#!/usr/bin/env python

'''
Add structural variants to a BAM.

The implementation lives in bamsurgeon.svengine.
'''

import argparse

from bamsurgeon.aligners import SUPPORTED_ALIGNERS
from bamsurgeon.svengine import main


def run():
    parser = argparse.ArgumentParser(description='adds SVs to reads, outputs modified reads as .bam along with mates')
    parser.add_argument('-v', '--varfile', dest='varFileName', required=True,
                        help='whitespace-delimited target regions for SV spike-in, see manual for syntax')
    parser.add_argument('-f', '--bamfile', dest='bamFileName', required=True,
                        help='sam/bam file from which to obtain reads')
    parser.add_argument('-r', '--reference', dest='refFasta', required=True,
                        help='reference genome, fasta indexed with bwa index _and_ samtools faidx')
    parser.add_argument('-o', '--outbam', dest='outBamFile', required=True,
                        help='.bam file name for output')
    parser.add_argument('-l', '--maxlibsize', dest='maxlibsize', default=600,
                        help="maximum fragment length of seq. library")
    parser.add_argument('-s', '--svfrac', dest='svfrac', default=1.0,
                        help="allele fraction of variant (default = 1.0)")
    parser.add_argument('--pad', default=None, type=int,
                        help='reference flank each side of a breakpoint (default = 2x --maxlibsize)')
    parser.add_argument('--readlen', default=None, type=int,
                        help='simulated read length (default: modal read length sampled from the BAM)')
    parser.add_argument('--mindepth', default=10, type=int,
                        help='minimum read depth in the breakend position to make mutation (default = 10)')
    parser.add_argument('--maxdepth', default=2000, type=int,
                        help='maximum read depth in the breakend position to make mutation (default = 2000)')
    parser.add_argument('--maxdfrac', default=0.1, type=float,
                        help='maximum discordant fraction (is_proper_pair / is_pair) of reads (default = 0.1)')
    parser.add_argument('--minctglen', dest='minctglen', default=None,
                        help=argparse.SUPPRESS)  # deprecated, see --pad
    parser.add_argument('-n', dest='maxmuts', default=None,
                        help="maximum number of mutations to make")
    parser.add_argument('-c', '--cnvfile', dest='cnvfile', default=None,
                        help="tabix-indexed list of genome-wide absolute copy number values (e.g. 2 alleles = no change)")
    parser.add_argument('--donorbam', dest='donorbam', default=None,
                        help='bam file supplying duplication interior reads; without it the interior is simulated')
    parser.add_argument('--ismean', dest='ismean', default=300,
                        help="mean insert size (default = estimate from region)")
    parser.add_argument('--issd', dest='issd', default=70,
                        help="insert size standard deviation (default = estimate from region)")
    parser.add_argument('--simerr', dest='simerr', default=0.0,
                        help='error rate for wgsim-generated reads')
    parser.add_argument('-p', '--procs', dest='procs', default=1,
                        help="split into multiple processes (default=1)")
    parser.add_argument('--inslib', default=None,
                        help='FASTA file containing library of possible insertions, use INS RND instead of INS filename to pick one')
    parser.add_argument('--aligner', default='backtrack',
                        help='supported aligners: ' + ','.join(sorted(SUPPORTED_ALIGNERS)))
    parser.add_argument('--alignopts', default=None,
                        help='aligner-specific options as comma delimited list of option1:value1,option2:value2,...')
    parser.add_argument('--alignerthreads', default=1,
                        help='threads used per realignment (default = 1)')
    parser.add_argument('--tagreads', action='store_true', default=False,
                        help='add BS tag to altered reads')
    parser.add_argument('--skipmerge', action='store_true', default=False,
                        help='do not merge spike-in reads back into original BAM')
    parser.add_argument('--keepsecondary', action='store_true', default=False,
                        help='keep secondary reads in final BAM')
    parser.add_argument('--debug', action='store_true', default=False,
                        help='output read tracking info to debug file, retain all intermediates')
    parser.add_argument('--maxopen', dest='maxopen', default=1000, type=int,
                        help="maximum number of open files during merge (default 1000)")
    parser.add_argument('--tmpdir', default='addsv.tmp',
                        help='temporary directory (default=addsv.tmp)')
    parser.add_argument('--seed', default=None, type=int,
                        help='seed random number generation')
    parser.add_argument('--allowN', action='store_true', default=False,
                        help='allow N in the reference slice, replace with A and warn (default: drop mutation)')
    parser.add_argument('--vcf', default='',
                        help="Path for the output VCF file. If not provided, the file will be saved in the current directory.")
    args = parser.parse_args()

    if args.pad is None:
        args.pad = 2 * int(args.maxlibsize)

    main(args)


if __name__ == '__main__':
    run()
