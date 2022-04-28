#!/usr/bin/env python

''' Adds a bogus read group to a readgroup-less BAM file '''

import os
import sys
import pysam
from uuid import uuid4
import logging

def modhead(header, rgid, fn):
    if 'RG' in header:
        logging.error("RG found in header, this script is not what you want!\n")
        sys.exit(usage())

    header['RG'] = [{'SM' : fn,
                     'LB' : 'bamsurgeon',
                     'CN' : 'BS',
                     'PU' : str(uuid4()),
                     'ID' : rgid,
                     'PL' : 'ILLUMINA' }]
    
    return header


def usage():
    return "usage: " + sys.argv[0] + " <BAM with no readgroups>"


if len(sys.argv) == 2:
    assert sys.argv[1].endswith('.bam'), usage()

    outbamfn = sys.argv[1].replace('.bam', '.BSRG.bam')
    rgid = str(uuid4())

    inbam  = pysam.AlignmentFile(sys.argv[1])
    outbam = pysam.AlignmentFile(outbamfn, 'wb', header=modhead(inbam.header, rgid, os.path.basename(outbamfn)))

    for read in inbam.fetch(until_eof=True):
        read.tags = read.tags + [('RG', rgid)]
        outbam.write(read)

    inbam.close()
    outbam.close()

else:
    sys.exit(usage())
