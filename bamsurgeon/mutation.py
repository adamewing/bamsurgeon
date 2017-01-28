#!/usr/bin/env python

from common import *
from collections import OrderedDict as od
import subprocess


def countBaseAtPos(bamfile,chrom,pos,mutid='null'):
    """ return list of bases at position chrom,pos
    """
    locstr = chrom + ":" + str(pos) + "-" + str(pos)
    args = ['samtools','mpileup',bamfile,'-r',locstr]

    p = subprocess.Popen(args,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    p.wait()
    pout = p.stdout.readlines()

    pileup = None 

    for line in pout:
        try:
            c = line.strip().split()
            assert len(c) > 5
            pileup = c[4].upper()
        except AssertionError:
            sys.stderr.write("INFO\t" + now() + "\t" + mutid + "\tmpileup failed, no coverage for base: " + chrom + ":" + str(pos) + "\n")
            return []
    bases = []
    if pileup:
        for b in pileup:
            if b in ['A','T','C','G']:
                bases.append(b)

    return bases


def makeins(read, start, ins, debug=False):
    assert len(read.seq) > len(ins) + 2
    
    if debug:
        print "DEBUG: INS: read.pos:", read.pos
        print "DEBUG: INS: start:   ", start
        print "DEBUG: INS: ins:     ", ins
        print "DEBUG: DEL: cigar:     ", read.cigarstring

    orig_len = len(read.seq)
    pos_in_read = start - read.pos + read.qstart


    if pos_in_read > 0: # insertion start in read span
        if debug:
            print "DEBUG: INS: pos_in_read:", pos_in_read
        left  = read.seq[:pos_in_read]
        right = read.seq[pos_in_read:]

        newseq = left + ins + right
        newseq = newseq[:orig_len]

    else: # insertion continues to the left of read
        right = read.seq[pos_in_read:]
        newseq = ins + right
        newseq = newseq[-orig_len:]

    if debug:
        print "DEBUG: INS: orig seq:", read.seq
        print "DEBUG: INS: newseq:  ", newseq
    return newseq


def makedel(read, chrom, start, end, ref, debug=False):
    assert len(read.seq) > end-start-2
    
    if debug:
        print "DEBUG: DEL: read.pos:", read.pos
        print "DEBUG: DEL: start:   ", start
        print "DEBUG: DEL: ins:     ", end
        print "DEBUG: DEL: cigar:     ", read.cigarstring
        print "DEBUG: DEL: orig seq:     ", read.seq

    orig_len = len(read.seq)
    #orig_end = read.pos + orig_len
    start_in_read = start - read.pos + read.qstart
    end_in_read = end - read.pos + read.qstart

    if debug:
        print "DEBUG: DEL: start_in_read:", start_in_read
        print "DEBUG: DEL: end_in_read:  ", end_in_read

    if start_in_read < 0: # deletion begins to the left of the read
        if debug:
            print "DEBUG: DEL: del begins to left of read." 

        assert end_in_read < orig_len
        right = read.seq[end_in_read:]
        left  = ref.fetch(chrom, start-(len(read.seq) - len(right)), start)

    elif end_in_read > orig_len: # deletion ends to the right of the read
        if debug:
            print "DEBUG: DEL: del ends to right of read."

        assert start_in_read > 0
        left  = read.seq[:start_in_read]
        right = ref.fetch(chrom, end, end+(len(read.seq) - len(left)))

    else:
        if debug:
            print "DEBUG: DEL: del starts and ends within read." 

        assert end_in_read <= orig_len and start_in_read >= 0 # deletion contained within the read
        left  = read.seq[:start_in_read]
        right = read.seq[end_in_read:]
        right += ref.fetch(chrom, read.pos+len(read.seq), read.pos+len(read.seq)+(len(read.seq)-len(left)-len(right)))

    if debug:
        print "DEBUG: DEL:  newseq:     ", left + right
    return left + right


def find_mate(read, bam):
    ''' AlignmentFile.mate() can return a non-primary alignment, so use this function instead '''
    chrom = read.next_reference_name
    for rec in bam.fetch(chrom, read.next_reference_start, read.next_reference_start+1):
        if rec.query_name == read.query_name and rec.reference_start == read.next_reference_start:
            if not rec.is_secondary and bin(rec.flag & 2048) != bin(2048):
                if rec.is_read1 != read.is_read1:
                    return rec
    return None


def mutate(args, log, bamfile, bammate, chrom, mutstart, mutend, mutpos_list, avoid=None, mutid_list=None, is_snv=False, mutbase_list=None, is_insertion=False, is_deletion=False, ins_seq=None, reffile=None, indel_start=None, indel_end=None):
    assert mutend > mutstart, "mutation start must occur before mutation end: " + mutid

    hasSNP = False

    outreads = od()
    mutreads = od()
    mutmates = od()

    region = 'haplo_' + chrom + '_' + str(mutstart) + '_' + str(mutend)

    maxfrac = None

    for pcol in bamfile.pileup(reference=chrom, start=mutstart, end=mutend, max_depth=int(args.maxdepth)):
        if pcol.pos:
            refbase = reffile.fetch(chrom, pcol.pos-1, pcol.pos)
            basepile = ''
            for pread in pcol.pileups:
                if avoid is not None and pread.alignment.qname in avoid:
                    print "WARN\t" + now() + "\t" + region + "\tdropped mutation due to read in --avoidlist", pread.alignment.qname
                    return True, False, {}, {}, {}

                # only consider primary alignments
                if pread.query_position is not None and not pread.alignment.is_secondary and bin(pread.alignment.flag & 2048) != bin(2048):
                    basepile += pread.alignment.seq[pread.query_position-1]
                    pairname = 'F' # read is first in pair
                    if pread.alignment.is_read2:
                        pairname = 'S' # read is second in pair
                    if not pread.alignment.is_paired:
                        pairname = 'U' # read is unpaired

                    extqname = ','.join((pread.alignment.qname,str(pread.alignment.pos),pairname))

                    if pcol.pos in mutpos_list:
                        if not pread.alignment.is_secondary and bin(pread.alignment.flag & 2048) != bin(2048) and not pread.alignment.mate_is_unmapped:
                            outreads[extqname] = pread.alignment
                            mutid = mutid_list[mutpos_list.index(pcol.pos)]

                            if is_snv:
                                if extqname not in mutreads:
                                    mutreads[extqname] = pread.alignment.seq

                                mutbase = mutbase_list[mutpos_list.index(pcol.pos)]
                                mutbases = list(mutreads[extqname])
                                mutbases[pread.query_position-1] = mutbase
                                mutread = ''.join(mutbases)
                                mutreads[extqname] = mutread

                            if is_insertion:
                                mutreads[extqname] = makeins(pread.alignment, indel_start, ins_seq)

                            if is_deletion:
                                mutreads[extqname] = makedel(pread.alignment, chrom, indel_start, indel_end, reffile)

                            mate = None
                            if not args.single:
                                try:
                                    mate = find_mate(pread.alignment, bammate)
                                except ValueError, e:
                                    raise ValueError('cannot find mate reference chrom for read %s, is this a single-ended BAM?' % pread.alignment.qname)

                                if mate is None:
                                    print "WARN\t" + now() + "\t" + mutid + "\twarning: no mate for", pread.alignment.qname
                                    if args.requirepaired:
                                        print "WARN\t" + now() + "\t" + mutid + "\tskipped mutation due to --requirepaired"
                                        return True, False, {}, {}, {}

                            if extqname not in mutmates:
                                mutmates[extqname] = mate

                            log.write(" ".join(('read',extqname,mutreads[extqname],"\n")))

                        if len(mutreads) > int(args.maxdepth):
                            sys.stderr.write("WARN\t" + now() + "\t" + mutid + "\tdepth at site is greater than cutoff, aborting mutation.\n")
                            return True, False, False, {}, {}, {}

            # make sure region doesn't have any changes that are likely SNPs
            # (trying to avoid messing with haplotypes)
            maxfrac = 0.0
            hasSNP  = False

            basepile = countBaseAtPos(args.bamFileName,chrom,pcol.pos,mutid=region)
            if basepile:
                majb = majorbase(basepile)
                minb = minorbase(basepile)

                frac = float(minb[1])/(float(majb[1])+float(minb[1]))
                if minb[0] == majb[0]:
                    frac = 0.0
                if frac > maxfrac:
                    maxfrac = frac
                if frac > float(args.snvfrac):
                    sys.stderr.write("WARN\t" + now() + "\t" + region + "\tdropped for proximity to SNP, nearby SNP MAF: " + str(frac)  + " (max snv frac: " + args.snvfrac + ")\n")
                    hasSNP = True
            else:
                sys.stderr.write("WARN\t" + now() + "\t" + region + "\tcould not pileup for region: " + chrom + ":" + str(pcol.pos) + "\n")
                if not args.ignorepileup:
                    hasSNP = True

    assert maxfrac is not None, "Error: could not pile up over region: %s" % region

    return False, hasSNP, maxfrac, outreads, mutreads, mutmates # todo: convert to class
