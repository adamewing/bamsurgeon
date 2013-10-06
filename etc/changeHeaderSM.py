#!/usr/bin/env python
 
import pysam
import sys
from re import sub
from os.path import basename
 
if len(sys.argv) == 2:
    assert sys.argv[1].endswith('.bam')
    newSM = sub('.bam$', '', basename(sys.argv[1]))
    outbamfn = sub('.bam$', '.replaceSM.bam', sys.argv[1])
    print "new SM:", newSM
    print "output BAM:", outbamfn
 
    inbam  = pysam.Samfile(sys.argv[1], 'rb')
 
    # change SM for all RG entries in header
    inheader = inbam.header
    for RG in inheader['RG']:
        RG['SM'] = newSM
 
    # create new BAM with modified header
    outbam = pysam.Samfile(outbamfn, 'wb', header=inheader)
 
    # put reads in new BAM
    n = 0
    for read in inbam.fetch(until_eof=True):
        outbam.write(read)
        n += 1
        if n % 10000000 == 0:
            print n,"reads replaced."
 
    inbam.close()
    outbam.close()
 
else:
    print "repaces 'SM' field in @RG section of header with name based on filename"
    print "usage:", sys.argv[0], "<.bam>"
