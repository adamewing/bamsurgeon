#!/usr/bin/env python

import sys
import pysam

import logging
FORMAT = '%(levelname)s %(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


if len(sys.argv) == 2:
    assert sys.argv[1].endswith('.bam')
    outfn = '.'.join(sys.argv[1].split('.')[:-1]) + '.BSonly.bam'

    logger.info('outputting reads with BS flag to %s' % outfn)

    inbam  = pysam.AlignmentFile(sys.argv[1], 'rb')
    outbam = pysam.AlignmentFile(outfn, 'wb', template=inbam)


    for read in inbam.fetch(until_eof=True):
        if read.has_tag('BS'):
            outbam.write(read)

    inbam.close()
    outbam.close()

else:
    sys.exit('usage: %s <bamsurgeon .bam file using --tagreads>' % sys.argv[0])
