#!/usr/bin/env python

import os
import pysam
import argparse
import logging

FORMAT = '%(levelname)s %(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

def fetchregions(infn, outfn, invcf, fasta_ref, window=1000):
    inbam  = pysam.AlignmentFile(infn, reference_filename=fasta_ref)
    outbam = pysam.AlignmentFile(outfn, 'wb', template=inbam)

    vcfh = pysam.VariantFile(invcf)

    for rec in vcfh:
        start = rec.pos - window
        # htslib folds INFO/END into rec.stop, so this covers symbolic SVs
        # too; the previous version subscripted the INFO value, which is a
        # scalar for Number=1 and would have raised
        end = max(rec.pos, rec.stop) + window

        if rec.chrom not in inbam.references:
            logger.warning("%s contig or chromosome not in %s" % (rec.chrom, infn))
            continue

        for read in inbam.fetch(rec.chrom, max(0, start), end):
            outbam.write(read)

    vcfh.close()

    inbam.close()
    outbam.close()


def main(args):
    assert args.bam.endswith('.bam'), "not a BAM file based on extension: " + args.bam
    assert os.path.exists(args.bam), "BAM file not found: " + args.bam
    assert os.path.exists(args.bam + '.bai'), "BAM index not found: " + args.bam + '.bai'
    assert os.path.exists(args.vcf), "VCF file not found: " + args.vcf
    assert args.vcf.endswith('.vcf') or args.vcf.endswith('.vcf.gz'), "not a VCF file based on extension " + args.vcf

    fetchregions(args.bam, args.out, args.vcf, args.fasta_ref, window=int(args.window))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="grab regions from BAM based on windows around records in a VCF file")
    parser.add_argument('-b', '--bam', required=True, help="BAM input (indexed)")
    parser.add_argument('-f', '--fasta-ref', help="FASTA reference (for CRAM)")
    parser.add_argument('-v', '--vcf', required=True, help="VCF input")
    parser.add_argument('-o', '--out', required=True, help="BAM output")
    parser.add_argument('-w', '--window', default=1000, help="window +/- VCF entry (default 1000)")
    args = parser.parse_args()
    main(args)
