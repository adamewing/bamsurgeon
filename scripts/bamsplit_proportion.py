#!/usr/bin/env python

import pysam
import sys
from re import sub
from random import random

import logging
FORMAT = '%(levelname)s %(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

if len(sys.argv) == 5:
    assert sys.argv[1].endswith('.bam')
    inbamfn = sys.argv[1]
    outbam1fn = sys.argv[3]
    outbam2fn = sys.argv[4]

    inbam = pysam.Samfile(inbamfn, 'rb')
    outbam1 = pysam.Samfile(outbam1fn, 'wb', template=inbam)
    outbam2 = pysam.Samfile(outbam2fn, 'wb', template=inbam)

    lastname = None
    lastread = None
    paired = False
    for read in inbam.fetch(until_eof=True):
        if not read.is_secondary:
            if read.qname == lastname:
                paired=True

            if paired:
                rnd = random()
                if rnd < float(sys.argv[2]):
                    outbam1.write(read)
                    outbam1.write(lastread)
                else:
                    outbam2.write(read)
                    outbam2.write(lastread)
                lastname = None
                lastread = None
                paired = False

            else:
                lastname = read.qname
                lastread = read
        else:
            logger.info("skipped secondary alignment: " + read.qname)

    outbam1.close()
    outbam2.close()
    inbam.close()
else:
    sys.exit("usage: %s <bam sorted by readname> <proportion split> <outbam1> <outbam2>" % sys.argv[0])

