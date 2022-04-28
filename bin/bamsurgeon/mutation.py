#!/usr/bin/env python

from bamsurgeon.common import *
from collections import OrderedDict as od
import subprocess

import logging
FORMAT = '%(levelname)s %(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def countBaseAtPos(bamfile,chrom,pos,mutid='null'):
    """ return list of bases at position chrom,pos
    """
    locstr = chrom + ":" + str(pos) + "-" + str(pos)
    args = ['samtools', 'mpileup', bamfile,'-Q', '0', '-r', locstr]

    p = subprocess.Popen(args,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    p.wait()
    pout = p.stdout.readlines()

    pileup = None 

    for line in pout:
        line = line.decode()
        try:
            c = line.strip().split()
            assert len(c) > 5
            pileup = c[4].upper()
        except AssertionError:
            logger.info(" mpileup failed, no coverage for base: " + chrom + ":" + str(pos))
            return []
    bases = []
    if pileup:
        for b in pileup:
            if b in ['A','T','C','G']:
                bases.append(b)

    return bases


def makeins(read, start, ins, debug=False):
    if len(read.seq) < len(ins) + 2:
        logger.warning("INDELs (ins) must be less than one read length, skipped read: %s" % read.query_name)
        return read.seq
    
    logger.debug("DEBUG: INS: read.pos: %d" % read.pos)
    logger.debug("DEBUG: INS: start: %d" % start)
    logger.debug("DEBUG: INS: ins: %s" % ins)
    logger.debug("DEBUG: INS: cigar: %s" % read.cigarstring)
    logger.debug("DEBUG: is_reverse: %s" % read.is_reverse)

    orig_len = len(read.seq)
    pos_in_read = None

    for (qpos, rpos) in read.get_aligned_pairs():
        if rpos == start:
            pos_in_read = qpos

    if pos_in_read is None:
        logger.warning("ref position %d not covered in read %s" % (start, read.query_name))
        return read.seq

    newseq = read.seq

    if pos_in_read > 0: # insertion start in read span
        logger.debug("DEBUG: INS: pos_in_read: %d" % pos_in_read)

        if not read.is_reverse:
            left  = read.seq[:pos_in_read]
            right = read.seq[pos_in_read:]

            newseq = left + ins + right
            newseq = newseq[:orig_len]

        else:
            pos_in_read = len(read.seq) - pos_in_read
            rcseq = rc(read.seq)

            left  = rcseq[:pos_in_read]
            right = rcseq[pos_in_read:]

            newseq = left + rc(ins) + right
            newseq = rc(newseq[:orig_len])

    logger.debug("DEBUG: INS: orig seq: %s" % read.seq)
    logger.debug("DEBUG: INS: newseq: %s" % newseq)

    return newseq


def makedel(read, chrom, start, end, ref, debug=False):
    if len(read.seq) < end-start-2:
        logger.warning("INDELs (del) must be less than one read length, skipped read: %s" % read.query_name)
        return read.seq
    
    logger.debug("DEBUG: DEL: read.pos: %d" % read.pos)
    logger.debug("DEBUG: DEL: start: %d" % start)
    logger.debug("DEBUG: DEL: ins: %d" % end)
    logger.debug("DEBUG: DEL: cigar: %s" % read.cigarstring)
    logger.debug("DEBUG: DEL: orig seq: %s" % read.seq)

    orig_len = len(read.seq)
    
    start_in_read = None
    end_in_read = None

    for (qpos, rpos) in read.get_aligned_pairs():
        if rpos == start:
            start_in_read = qpos

        if rpos == end:
            end_in_read = qpos

    if start_in_read is None and read.get_reference_positions()[0] > start:
        start_in_read = start-read.get_reference_positions()[0]

    if end_in_read is None and read.get_reference_positions()[-1] < end:
        end_in_read = orig_len + (end-read.get_reference_positions()[-1])

    if start_in_read is None:
        logger.warning("ref position %d not covered in read %s" % (start, read.query_name))
        return read.seq

    if end_in_read is None:
        logger.warning("ref position %d not covered in read %s" % (end, read.query_name))
        return read.seq

    logger.debug("DEBUG: DEL: start_in_read: %d" % start_in_read)
    logger.debug("DEBUG: DEL: end_in_read: %d" % end_in_read)

    if start_in_read < 0: # deletion begins to the left of the read
        
        logger.debug("DEBUG: DEL: del begins to left of read.")

        assert end_in_read < orig_len
        right = read.seq[end_in_read:]
        left  = ref.fetch(chrom, start-(len(read.seq) - len(right)), start)

    elif end_in_read > orig_len: # deletion ends to the right of the read
        logger.debug("DEBUG: DEL: del ends to right of read.")

        assert start_in_read > 0
        left  = read.seq[:start_in_read]
        right = ref.fetch(chrom, end, end+(len(read.seq) - len(left)))

    else:
        logger.debug("DEBUG: DEL: del starts and ends within read.") 

        assert end_in_read <= orig_len and start_in_read >= 0 # deletion contained within the read
        left  = read.seq[:start_in_read]
        right = read.seq[end_in_read:]
        right += ref.fetch(chrom, read.pos+len(read.seq), read.pos+len(read.seq)+(len(read.seq)-len(left)-len(right)))

    if debug:
        logger.debug("DEBUG: DEL: newseq: %s" % (left + right))

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
    assert mutend > mutstart, "mutation start must occur before mutation end"

    hasSNP = False

    outreads = od()
    mutreads = od()
    mutmates = od()

    region = 'haplo_' + chrom + '_' + str(mutstart) + '_' + str(mutend)

    maxfrac = None

    for pcol in bamfile.pileup(reference=chrom, start=mutstart-1, end=mutend+1, max_depth=int(args.maxdepth), ignore_overlaps=False):
        if pcol.pos:
            if args.ignorepileup and (pcol.pos < mutstart-1 or pcol.pos > mutend+1):
                continue

            refbase = reffile.fetch(chrom, pcol.pos-1, pcol.pos)
            basepile = ''
            for pread in pcol.pileups:
                if avoid is not None and pread.alignment.qname in avoid:
                    logger.warning(region + " dropped mutation due to read in --avoidlist " + pread.alignment.qname)
                    return True, False, maxfrac, {}, {}, {}

                # only consider primary alignments
                if pread.query_position is not None and not pread.alignment.is_secondary and bin(pread.alignment.flag & 2048) != bin(2048):
                    basepile += pread.alignment.seq[pread.query_position-1]
                    pairname = 'F' # read is first in pair
                    if pread.alignment.is_read2:
                        pairname = 'S' # read is second in pair
                    if not pread.alignment.is_paired:
                        pairname = 'U' # read is unpaired

                    extqname = ','.join((pread.alignment.qname,str(pread.alignment.pos),pairname))

                    if pcol.pos+1 in mutpos_list:

                        if not pread.alignment.is_secondary and bin(pread.alignment.flag & 2048) != bin(2048) and not pread.alignment.mate_is_unmapped:
                            outreads[extqname] = pread.alignment
                            mutid = mutid_list[mutpos_list.index(pcol.pos+1)]

                            if is_snv:
                                if extqname not in mutreads:
                                    mutreads[extqname] = pread.alignment.seq

                                mutbase = mutbase_list[mutpos_list.index(pcol.pos+1)]
                                mutbases = list(mutreads[extqname])
                                mutbases[pread.query_position] = mutbase
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
                                except ValueError:
                                    raise ValueError('cannot find mate reference chrom for read %s, is this a single-ended BAM?' % pread.alignment.qname)

                                if mate is None:
                                    logger.warning(mutid + " warning: no mate for " + pread.alignment.qname)
                                    if args.requirepaired:
                                        logger.warning(mutid + " skipped mutation due to --requirepaired")
                                        return True, False, {}, {}, {}

                            if extqname not in mutmates:
                                mutmates[extqname] = mate

                            log.write(" ".join(('read',extqname,mutreads[extqname],"\n")))

                        if len(mutreads) > int(args.maxdepth):
                            logger.warning("depth at site is greater than cutoff, aborting mutation")
                            return True, False, maxfrac, {}, {}, {}

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
                    logger.warning(region + " dropped for proximity to SNP, nearby SNP MAF: " + str(frac)  + " (max snv frac: " + args.snvfrac + ")")
                    hasSNP = True
            else:
                logger.warning(region + " could not pileup for region: " + chrom + ":" + str(pcol.pos))
                if not args.ignorepileup:
                    hasSNP = True

    if maxfrac is None:
        logger.warning("could not pile up over region: %s" % region)
        return True, False, maxfrac, {}, {}, {}

    return False, hasSNP, maxfrac, outreads, mutreads, mutmates # todo: convert to class
