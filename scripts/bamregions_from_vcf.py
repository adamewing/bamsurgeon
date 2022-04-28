#!/usr/bin/env python

import vcf
import os
import sys
import pysam
import argparse
import logging

def fetchregions(infn, outfn, invcf, window=1000):
    inbam  = pysam.AlignmentFile(infn)
    outbam = pysam.AlignmentFile(outfn, 'wb', template=inbam)

    vcfh = vcf.Reader(filename=invcf)

    for rec in vcfh:
        start = rec.POS - window
        end   = rec.POS + window

        if 'END' in rec.INFO:
            end = int(rec.INFO.get('END')[0]) + window

        if rec.CHROM not in inbam.references:
            logging.warn("WARNING: " + rec.CHROM + " contig or chromosome not in " + infn + "\n")
            continue

        for read in inbam.fetch(rec.CHROM, start, end):
            outbam.write(read)

    inbam.close()
    outbam.close()


def main(args):
    assert args.bam.endswith('.bam'), "not a BAM file based on extension: " + args.bam
    assert os.path.exists(args.bam), "BAM file not found: " + args.bam
    assert os.path.exists(args.bam + '.bai'), "BAM index not found: " + args.bam + '.bai'
    assert os.path.exists(args.vcf), "VCF file not found: " + args.vcf
    assert args.vcf.endswith('.vcf') or args.vcf.endswith('.vcf.gz'), "not a VCF file based on extension " + args.vcf

    fetchregions(args.bam, args.out, args.vcf, window=int(args.window))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="grab regions from BAM based on windows around records in a VCF file")
    parser.add_argument('-b', '--bam', required=True, help="BAM input (indexed)")
    parser.add_argument('-v', '--vcf', required=True, help="VCF input")
    parser.add_argument('-o', '--out', required=True, help="BAM output")
    parser.add_argument('-w', '--window', default=1000, help="window +/- VCF entry (default 1000)")
    args = parser.parse_args()
    main(args)

