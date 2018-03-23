#!/usr/bin/env python

import argparse
import pysam
import sys
import datetime

from subprocess import call
from re import sub
from os.path import basename
from os import remove, rename
from uuid import uuid4


def now():
    return str(datetime.datetime.now())


def getRG(tags):
    ''' fetch RG, tags is from a pysam.AlignedRead.tags, returns RG name '''
    for tag, val in tags:
        if tag == 'RG':
            return val
    return None


def putRG(tags, rg):
    ''' replace RG, tags is from a pysam.AlignedRead.tags, returns list of tags '''
    out = []
    for tag, val in tags:
        if tag == 'RG':
            out.append((tag, rg))
        else:
            out.append((tag, val))
    return out


def samrec(read, bam, IDRG, newname=None):
    ''' output sam formatted record from pysam.AlignedRead '''
    fields = []

    if newname is not None:
        fields.append(newname)
    else:
        fields.append(read.qname)

    fields.append(str(read.flag))

    if read.is_unmapped:
        if read.mate_is_unmapped:
            fields.append('*')
        else:
            fields.append(bam.getrname(read.rnext))
    else:
        fields.append(bam.getrname(read.tid))

    if read.is_unmapped:
        if read.mate_is_unmapped:
            fields.append('0') # was 0
        else:
            fields.append(str(read.mpos+1)) #avoid OBOE, pysam is always 0-based, SAM is 1-based, BAM is 0-based.
        fields.append('0') # was 0
    else:
        fields.append(str(read.pos+1)) #avoid OBOE
        fields.append(str(read.mapq))

    # unmapped reads should have '*' as CIGAR string
    if read.is_unmapped:
        fields.append('*')
    else:
        fields.append(read.cigarstring)

    if read.tid == read.rnext:
        fields.append('=')
    else:
        if read.is_unmapped or read.mate_is_unmapped:
            fields.append('=')
        else:
            fields.append(bam.getrname(read.rnext))

    if read.mate_is_unmapped:
        if read.is_unmapped:
            fields.append('0') # was 0
        else:
            fields.append(str(read.pos+1)) # avoid OBOE
    else:
        fields.append(str(read.mpos+1)) # avoid OBOE

    if read.is_unmapped or read.mate_is_unmapped:
        fields.append('0')
    else:
        fields.append(str(read.isize))
    fields.append(read.seq)
    fields.append(read.qual)
    
    for tag in read.tags:
        tagname, tagval = tag
        # retain only certain tags
        if tagname in ('NM', 'MD', 'AS', 'XS', 'RG'):
            tagtype = None

            # note, python has no 'char' type so tags with 'A' type will be converted to 'Z'
            if type(tagval) == type('a'):
                tagtype = 'Z'
            if type(tagval) == type(1):
                tagtype = 'i'
            if type(tagval) == type(1.0):
                tagtype = 'f'

            if tagname == 'RG':
                if tagval not in IDRG:
                    sys.stderr.write("ERROR:\t" + now() + "\tread group: " + tagval + " not found!\n")
                assert tagval in IDRG # make sure readgroup is value (i.e. in the header)
                tagval = IDRG[tagval]

            if tagtype is not None:
                fields.append(':'.join((tagname, tagtype, str(tagval))))

    return '\t'.join(fields)


def makebam(sam, fai, threads, mem):
    ''' sam --> sorted bam '''
    outbam = sub('.sam$', '.bam', sam)
    cmd = ['samtools', 'view', '-bt', fai, '-o', outbam, sam]
    sys.stderr.write(sam + ' --> ' + outbam + ': ' + ' '.join(cmd) + '\n')
    call(cmd)
    #remove(sam)

    outsort = sub('.bam$', '.sorted.bam', outbam)
    cmd = ['samtools', 'sort', '-m', str(mem), '-@', str(threads), '-T', outsort, '-o',  outsort, outbam]
    outsort += '.bam'
    sys.stderr.write(outbam + ' --> ' + outsort + ': ' + ' '.join(cmd) + '\n')
    call(cmd)

    #remove(outbam)


def main(args):
    assert args.bam[0].endswith('.bam')
    assert args.fai.endswith('.fai')
    outsamfn = sub('.bam$', '.postprocessed.sam', args.bam[0])

    bam = pysam.Samfile(args.bam[0], 'rb')

    PURG = {}
    IDRG = {}
    header = bam.header

    if 'RG' in header:
        newSM = sub('.bam$', '', basename(args.bam[0]))
        for RG in header['RG']:
            RG['SM'] = newSM
            RG['LB'] = 'bamsurgeon'
            RG['CN'] = 'BS'
            if 'PU' in RG and RG['PU'] not in PURG:
                PU = str(uuid4())
                PURG[RG['PU']] = PU
                RG['PU'] = PU
            if 'ID' in RG and RG['ID'] not in IDRG:
                ID = str(uuid4())
                IDRG[RG['ID']] = ID
                RG['ID'] = ID 

    if 'PG' in header:
        del header['PG']
        header['PG'] = [{'ID': 'bamsurgeon', 'PN': 'bamsurgeon'}]

    outsam = pysam.Samfile(outsamfn, 'wh', header=header)
    outsam.close()

    paired = {} # track read pairs

    # counters for debug
    n = 0 # number of reads
    p = 0 # numbar of paired reads
    u = 0 # number of unpaired reads
    w = 0 # reads written
    m = 0 # mates found

    fixed_strand  = 0
    fixed_rg_pair = 0
    fixed_matepos = 0
    fixed_tlen    = 0
    fixed_unmap   = 0
    fixed_materef = 0

    tick = 100000
    try:
        tick = int((bam.mapped + bam.unmapped) * 0.01)
        if tick == 0:
            tick = 1
        sys.stderr.write("INFO\t" + now() + "\toutputting status every " + str(tick) + " reads (1%) ...\n")
    except ValueError as e:
        sys.stderr.write("INFO\t" + now() + "\tno index found, outputting status every " + str(tick) + " reads.\n")

    outsam = open(outsamfn, 'a')

    for read in bam.fetch(until_eof=True):
        n += 1
        if read.is_paired and not read.is_secondary:
            p += 1
            if read.qname in paired:
                # make sure paired read groups match
                rg = getRG(read.tags)
                if rg != getRG(paired[read.qname].tags):
                    read.tags = putRG(read.tags, rg)
                    paired[read.qname].tags = putRG(paired[read.qname].tags, rg)

                    assert rg == getRG(paired[read.qname].tags)
                    fixed_rg_pair += 1

                # fix strand 
                if read.mate_is_reverse != paired[read.qname].is_reverse or paired[read.qname].mate_is_reverse != read.is_reverse:
                    read.mate_is_reverse = paired[read.qname].is_reverse
                    paired[read.qname].mate_is_reverse = read.is_reverse

                    assert read.mate_is_reverse == paired[read.qname].is_reverse and paired[read.qname].mate_is_reverse == read.is_reverse
                    fixed_strand += 1

                # fix mate position
                if read.pnext != paired[read.qname].pos or paired[read.qname].pnext != read.pos:
                    read.pnext = paired[read.qname].pos
                    paired[read.qname].pnext = read.pos

                    assert read.pnext == paired[read.qname].pos and paired[read.qname].pnext == read.pos
                    fixed_matepos += 1

                # fix unmapped flag
                if read.mate_is_unmapped != paired[read.qname].is_unmapped or paired[read.qname].mate_is_unmapped != read.is_unmapped:
                    read.mate_is_unmapped = paired[read.qname].is_unmapped
                    paired[read.qname].mate_is_unmapped = read.is_unmapped

                    assert read.mate_is_unmapped == paired[read.qname].is_unmapped and paired[read.qname].mate_is_unmapped == read.is_unmapped
                    fixed_unmap += 1

                # fix mate ref
                if read.tid != paired[read.qname].rnext or paired[read.qname].tid != read.rnext:
                    read.rnext = paired[read.qname].tid
                    paired[read.qname].rnext = read.tid

                    assert read.tid == paired[read.qname].rnext and paired[read.qname].tid == read.rnext
                    fixed_materef += 1

                # fix tlen (left - (right + read length) where left < right)
                if not read.is_unmapped and not paired[read.qname].is_unmapped and read.tid == paired[read.qname].tid:
                    if abs(read.tlen) != abs(paired[read.qname].tlen):
                        read.tlen = min(read.pos, read.pnext)-(max(read.pos, read.pnext)+read.rlen)
                        paired[read.qname].tlen = 0-read.tlen
                        
                        assert abs(read.tlen) == abs(paired[read.qname].tlen)
                        fixed_tlen += 1

                newname = None 
                if args.rename:
                    newname = str(uuid4())

                outsam.write(samrec(read, bam, IDRG, newname=newname) + '\n')               # output read
                outsam.write(samrec(paired[read.qname], bam, IDRG, newname=newname) + '\n') # output mate
                del paired[read.qname]
                w += 1
                m += 1
            else:
                paired[read.qname] = read
                w += 1
        else:
            if not read.is_secondary:
                u += 1

                newname = None
                if args.rename:
                    newname = str(uuid4())

                outsam.write(samrec(read, bam, IDRG, newname=newname) + '\n')
                w += 1

        if n % tick == 0:
            sys.stderr.write('\t'.join(map(str, ('processed',n,'reads:',p,'paired',u,'unpaired',w,'written',m,'mates found.'))) + '\n')
            sys.stderr.write('\t'.join(map(str, ('fixed strand:', fixed_strand, 'fixed RG pair:', fixed_rg_pair, 'fixed mate pos:', fixed_matepos))) + '\n')
            sys.stderr.write('\t'.join(map(str, ('fixed unmapped flag:', fixed_unmap, 'fixed mate ref:', fixed_materef, 'fixed tlen:', fixed_tlen))) + '\n')

    if len(paired.keys()) > 0:
        sys.stderr.write("WARNING:\t" + now() + "\tfound " + str(len(paired.keys())) + " orphaned paired reads that were not output!\n") 

    outsam.close()
    makebam(outsamfn, args.fai, args.threads, args.mem)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Postprocess BAM files generated by bamsurgeon to ensure consistency and compliance with SAM spec.')
    parser.add_argument(metavar='<BAM file>', dest='bam', nargs=1, help='BAM file (from bamsurgeon output)')
    parser.add_argument('-f', '--fai', dest='fai', required=True, help='.fai index, generated with samtools faidx on reference FASTA')
    parser.add_argument('-t', '--sort-threads', dest='threads', default=1, help='threads for sorting with samtools (-@)')
    parser.add_argument('-m', '--sort-mem', dest='mem', default='4G', help='memory PER THREAD for sorting with samtools (-m)')
    parser.add_argument('--rename', action='store_true', default=False, help='rename reads to uuids')
    args = parser.parse_args()
    main(args)
