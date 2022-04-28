#!/usr/bin/env python

import os
import pysam
import argparse
import logging
import subprocess

from collections import OrderedDict as od

logger = logging.getLogger(__name__)
FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger.setLevel(logging.INFO)


def main(args):
    assert os.path.exists(args.fasta + '.fai'), 'please run samtools faidx %s' % args.fasta

    fa  = pysam.FastaFile(args.fasta)
    bam = pysam.AlignmentFile(args.bam)

    bam_reflen = od([(ref, length) for ref, length in zip(bam.references, bam.lengths)])

    with open(args.outfa, 'w') as outfa:
        for ref in bam.references:
            assert ref in fa.references, 'reference not shared: %s' % ref
            assert fa.get_reference_length(ref) == bam_reflen[ref], 'length mismatch on contig/chr %s' % ref

            logger.info('writing %s to %s (len = %d)' % (ref, args.outfa, bam_reflen[ref]))

            outfa.write('>%s\n%s\n' % (ref, fa.fetch(ref)))

    


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="re-arrange indexed FASTA to match BAM reference order")
    parser.add_argument('-f', '--fasta', required=True, help="FASTA must be indexed with samtools faidx")
    parser.add_argument('-b', '--bam', required=True, help="BAM")
    parser.add_argument('-o', '--outfa', default='bam_matched.fa', help='output FASTA (default=bam_matched.fa)')
    args = parser.parse_args()
    main(args)
