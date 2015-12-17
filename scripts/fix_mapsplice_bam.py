#!/usr/bin/env python

''' Script for correcting broken mapsplice BAMs using pysam
    Adam Ewing (ewingad@soe.ucsc.edu)
'''

import sys
import pysam
import os
import subprocess
from re import sub

def namesort_bam(bamfile):
    sortbase = bamfile + ".namesort"
    sortfn   = sortbase + ".bam"
    sortargs = ['samtools','sort','-n','-@','8','-m','2G',bamfile,'-o',sortbase]
    print "sorting, cmd: " + " ".join(sortargs)
    subprocess.call(sortargs)
    return sortfn

def possort_bam(bamfile):
    sortbase = bamfile + ".sort"
    sortfn   = sortbase + ".bam"
    sortargs = ['samtools','sort','-@','8','-m','2G',bamfile,'-o',sortbase]
    print "sorting, cmd: " + " ".join(sortargs)
    subprocess.call(sortargs)
    os.rename(sortfn,bamfile)

def validate(reads):
    if len(reads) > 2:
        return False
    if len(reads) == 2:
        if reads[0].is_unmapped and reads[1].is_unmapped:
            return True

        has_read1 = has_read2 = False

        for read in reads:
            if read.is_read1 and read.is_read2:
                return False
            if not read.mate_is_unmapped and read.mpos < 0:
                return False

            if read.is_read1:
                has_read1 = True
            if read.is_read2:
                has_read2 = True

        if not (has_read1 and has_read2):
            return False

        # check isize
        if abs(reads[0].isize) != abs(reads[1].isize):
            return False

        # check paired flag
        if not (reads[0].is_paired and reads[1].is_paired):
            return False

        # mate strand agreement
        if reads[0].is_reverse != reads[1].mate_is_reverse:
            return False

        if reads[1].is_reverse != reads[0].mate_is_reverse:
            return False

        # mate position and refid agreement
        if reads[0].tid != reads[1].rnext or reads[1].tid != reads[0].rnext:
            return False

        if reads[0].pos != reads[1].mpos or reads[1].pos != reads[0].mpos:
            return False

    return True

def fixmates(reads):
    ''' if there are more than 2 reads in a pair:
        1. if one is marked non-primary, remove it
        2. if (1) results in two 'read1' reads mark one as read
        3. TODO: try to base (2) on FR read mapping scheme '''

    if len(reads) > 2:
        newreads = []
        for read in reads:
            if not read.is_secondary:
                newreads.append(read)
        if len(newreads) == 2:
            if validate(newreads): # if rejecting a non-primary alignment fixes it
                return newreads
        else: # unsalvagable at present, don't output
            sys.stderr.write("rejecting (three primary alignments for pair):" + newreads[0].qname + "\n")
            return []
        reads = newreads

    if len(reads) == 2:
        # fix mate strand agreement, position, refid
        reads[1].mate_is_reverse = reads[0].is_reverse
        reads[0].mate_is_reverse = reads[1].is_reverse

        reads[0].rnext = reads[1].tid
        reads[1].rnext = reads[0].tid

        reads[0].mpos = reads[1].pos
        reads[1].mpos = reads[0].pos

        if validate(reads):
            return reads

        # fix unpaired by flag
        if not reads[0].is_paired:
            reads[0].is_paired=True
        if not reads[1].is_paired:
            reads[1].is_paired=True
        if validate(reads):
            return reads

        # fix one-end anchored
        if (not reads[0].is_unmapped) and reads[1].is_unmapped:
            reads[0].mate_is_unmapped = True
            if validate(reads):
                return reads

        elif (not reads[1].is_unmapped) and reads[0].is_unmapped:
            reads[1].mate_is_unmapped = True
            if validate(reads):
                return reads

        # fix insert size, set to smallest of the pair (mapsplice sets huge isizes sometimes)
        if abs(reads[0].isize) != abs(reads[1].isize):
            if (not reads[0].is_unmapped) and (not reads[1].is_unmapped) and (reads[0].tid == reads[1].tid):
                reads[0].isize = reads[0].pos - reads[1].pos
                reads[1].isize = reads[1].pos - reads[0].pos
            else:
                reads[0].isize = reads[1].isize = 0

        # fix if mates don't have position/tid set
        if (not reads[0].mate_is_unmapped and reads[0].mpos < 0) or (not reads[1].mate_is_unmapped and reads[1].mpos < 0):
            reads[0].mpos = reads[1].pos
            reads[0].rnext = reads[1].tid
            reads[1].mpos = reads[0].pos
            reads[1].rnext = reads[0].tid

        if validate(reads):
            return reads

        # try to fix both-ended reads (where a read is marked both end1 and end2)
        newreads = []
        for read in reads:
            if read.is_read1 and read.is_read2:
                # try to infer correct order from query name
                if read.qname.endswith('/1'):
                    read.is_read2 = False
                elif read.qname.endswith('/2'):
                    read.is_read1 = False
            newreads.append(read)
        if validate(newreads): # if fixing both-ended reads is enough
            return newreads
        else:
            reads = newreads

        # fix situations where there is no read1 or no read2 in a pair
        has_read1 = has_read2 = False
        for read in reads:
            if read.is_read1:
                has_read1 = True
            if read.is_read2:
                has_read2 = True

        # try to assign based on qname
        newreads = []
        if not (has_read1 and has_read2):
            for read in reads:
                if read.qname.endswith('/1'):
                    read.is_read1 = True
                    read.is_read2 = False
                elif read.qname.endswith('/2'):
                    read.is_read2 = True
                    read.is_read1 = False
                newreads.append(read)
            reads = newreads

        if validate(reads):
            return reads
        else:
            # arbitrary assignment as last ditch option
            reads[0].is_read1 = True
            reads[0].is_read2 = False
            reads[1].is_read1 = False
            reads[1].is_read2 = True
            if validate(reads):
                return reads

    sys.stderr.write("rejecting (could not correct):" + reads[0].qname + "\n")
    return []


def writereads(reads, outbam):
    for read in reads:
        read.qname = sub('/[12]', '', read.qname)
        outbam.write(read)

if len(sys.argv) == 2:
    assert sys.argv[1].endswith('bam')
    inbamfile  = namesort_bam(sys.argv[1])
    outbamfile = sub('bam$','fix.bam',sys.argv[1])

    inbam  = pysam.Samfile(inbamfile, 'rb')
    outbam = pysam.Samfile(outbamfile, 'wb', template=inbam)

    reads = []
    passed = fixed = 0
    lastname = None

    for read in inbam.fetch(until_eof=True):
        basename = sub('/[12]', '', read.qname)
        if basename == lastname:
            reads.append(read)
        else:
            if validate(reads):
                passed += 1
                writereads(reads, outbam)
            else:
                reads = fixmates(reads)
                writereads(reads, outbam)
                fixed += 1

            reads = []
            lastname = basename
            reads.append(read)

    # handle last reads
    if validate(reads):
        passed += 1
        writereads(reads, outbam)
    else:
        reads = fixmates(reads)
        writereads(reads, outbam)
        fixed += 1

    outbam.close()
    os.remove(inbamfile)
    print 'groups passed:',passed,'fixed:',fixed,'... sorting'
    possort_bam(outbamfile)
else:
    print "corrects mapsplice .bam files.\n"
    print "usage:", sys.argv[0], "<BAM generated by mapsplice>"

