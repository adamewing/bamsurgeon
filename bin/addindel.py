#!/usr/bin/env python

'''
Add small insertions and deletions to a BAM.

The implementation lives in bamsurgeon.snvindel, which addsnv.py shares.
'''

import os
import logging
import argparse

from bamsurgeon.aligners import SUPPORTED_ALIGNERS
from bamsurgeon.snvindel import read_indel_targets, run_spikein
from bamsurgeon.varinput import sample_read_length

logger = logging.getLogger(__name__)


def main(args):
    read_length = sample_read_length(args.bamFileName, args.refFasta) or 0
    clusters = read_indel_targets(args.varFileName, int(args.numsnvs), read_length)
    run_spikein(args, clusters, 'addindel')


def run():
    parser = argparse.ArgumentParser(description='adds INDELs to reads, outputs modified reads as .bam along with mates')
    parser.add_argument('-v', '--varfile', dest='varFileName', required=True, help='Target regions to try and add an indel, as BED')
    parser.add_argument('-f', '--bamfile', dest='bamFileName', required=True, help='sam/bam file from which to obtain reads')
    parser.add_argument('-r', '--reference', dest='refFasta', required=True, help='reference genome, fasta indexed with bwa index _and_ samtools faidx')
    parser.add_argument('-o', '--outbam', dest='outBamFile', required=True, help='.bam file name for output')
    parser.add_argument('-s', '--snvfrac', dest='snvfrac', default=1, help='maximum allowable linked SNP MAF (for avoiding haplotypes) (default = 1)')
    parser.add_argument('-m', '--mutfrac', dest='mutfrac', default=0.5, help='allelic fraction at which to make indels (default = 0.5)')
    parser.add_argument('-n', '--numsnvs', dest='numsnvs', default=0, help="maximum number of mutations to try (default: entire input)")
    parser.add_argument('-c', '--cnvfile', dest='cnvfile', default=None, help="tabix-indexed list of genome-wide absolute copy number values (e.g. 2 alleles = no change)")
    parser.add_argument('-d', '--coverdiff', dest='coverdiff', default=0.1, help="allow difference in input and output coverage (default=0.1)")
    parser.add_argument('-p', '--procs', dest='procs', default=1, help="split into multiple processes (default=1)")
    parser.add_argument('--picardjar', default=os.environ.get('BAMSURGEON_PICARD_JAR'), help='path to picard.jar (default: $BAMSURGEON_PICARD_JAR)')
    parser.add_argument('--mindepth', default=10, help='minimum read depth to make mutation (default = 10)')
    parser.add_argument('--maxdepth', default=2000, help='maximum read depth to make mutation (default = 2000)')
    parser.add_argument('--minmutreads', default=3, help='minimum number of mutated reads to output per site')
    parser.add_argument('--avoidreads', default=None, help='file of read names to avoid (mutations will be skipped if overlap)')
    parser.add_argument('--nomut', action='store_true', default=False, help="dry run")
    parser.add_argument('--det', action='store_true', default=False, help=argparse.SUPPRESS)  # deprecated, no effect
    parser.add_argument('--ignoresnps', action='store_true', default=False, help="skip the scan for confounding SNPs sharing reads with the mutation")
    parser.add_argument('--insane', action='store_true', default=False, help="skip the coverage sanity check comparing input and output depth")
    parser.add_argument('--force', action='store_true', default=False, help="force mutation to happen regardless of nearby SNP or low coverage")
    parser.add_argument('--single', action='store_true', default=False, help="input BAM is simgle-ended (default is paired-end)")
    parser.add_argument('--maxopen', dest='maxopen', default=1000, help="maximum number of open files during merge (default 1000)")
    parser.add_argument('--requirepaired', action='store_true', default=False, help='skip mutations if unpaired reads are present')
    parser.add_argument('--aligner', default='backtrack', help='supported aligners: ' + ','.join(sorted(SUPPORTED_ALIGNERS)))
    parser.add_argument('--alignopts', default=None, help='aligner-specific options as comma delimited list of option1:value1,option2:value2,...')
    parser.add_argument('--alignerthreads', default=1, help='threads used per realignment (default = 1)')
    parser.add_argument('--tagreads', action='store_true', default=False, help='add BS tag to altered reads')
    parser.add_argument('--skipmerge', action='store_true', default=False, help="final output is tmp file to be merged")
    parser.add_argument('--ignorepileup', action='store_true', default=False, help="do not check pileup depth in mutation regions")
    parser.add_argument('--tmpdir', default='addindel.tmp', help='temporary directory (default=addindel.tmp)')
    parser.add_argument('--seed', default=None, type=int, help='seed random number generation')
    parser.add_argument('--vcf', default='', help="Path for the output VCF file. If not provided, the file will be saved in the current directory.")
    args = parser.parse_args()

    # addindel never took these, but the shared implementation reads them
    args.haplosize = 0
    args.ignoreref = False

    if args.det:
        logger.warning('--det is deprecated and has no effect')

    if args.picardjar is None:
        parser.error('--picardjar is required, or set BAMSURGEON_PICARD_JAR in the environment')

    main(args)


if __name__ == '__main__':
    run()
