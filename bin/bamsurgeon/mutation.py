#!/usr/bin/env python

from bamsurgeon.common import *
from collections import OrderedDict as od
from dataclasses import dataclass, field
import subprocess

import logging
FORMAT = '%(levelname)s %(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


# pysam's pileup() default. The pileups are now built without it so that read
# selection sees every read, and it is applied here instead, where a base call
# too poor to trust is actually the question being asked.
SNP_SCAN_MIN_BASEQUAL = 13


def column_bases(pcol, min_basequal=SNP_SCAN_MIN_BASEQUAL):
    """
    ACGT bases at a pileup column, from primary alignments only.

    Replaces countBaseAtPos(), which forked `samtools mpileup` once per
    column -- hundreds of process spawns per mutation, and the dominant cost
    of addsnv/addindel. The data is already in the column being iterated.

    That function also passed pysam's 0-based pcol.pos into a 1-based mpileup
    region string, so it reported the bases one position to the left of the
    column it was called for.
    """
    bases = []

    for pread in pcol.pileups:
        if pread.is_del or pread.is_refskip or pread.query_position is None:
            continue

        aln = pread.alignment
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary \
                or aln.is_duplicate or aln.is_qcfail:
            continue

        quals = aln.query_qualities
        if quals is not None and quals[pread.query_position] < min_basequal:
            continue

        base = aln.query_sequence[pread.query_position].upper()
        if base in ('A', 'C', 'G', 'T'):
            bases.append(base)

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


@dataclass
class MutateResult:
    """
    Outcome of collecting and editing the reads over a mutation region.

    A dataclass rather than a tuple because one of the six former return
    sites returned five values, so `--requirepaired` with a mateless read
    unpacked short and raised instead of skipping the site.
    """
    failed: bool = False
    has_snp: bool = False
    maxfrac: float = None
    outreads: dict = field(default_factory=od)
    mutreads: dict = field(default_factory=od)
    mutmates: dict = field(default_factory=od)


def mutate(args, log, bamfile, bammate, chrom, mutstart, mutend, mutpos_list, avoid=None, mutid_list=None, is_snv=False, mutbase_list=None, is_insertion=False, is_deletion=False, ins_seq=None, reffile=None, indel_start=None, indel_end=None):
    assert mutend > mutstart, "mutation start must occur before mutation end"

    result = MutateResult()
    region = 'haplo_' + chrom + '_' + str(mutstart) + '_' + str(mutend)

    snvfrac = float(args.snvfrac)

    # --ignoresnps skips the scan for confounding nearby SNPs entirely, rather
    # than running it and discarding the answer
    ignoresnps = getattr(args, 'ignoresnps', False)

    # min_base_quality=0 because this pileup decides which reads get edited,
    # and how confident the sequencer was in one base call has no bearing on
    # that. pysam's default of 13 hid 22 of the 29 reads over a site in the
    # ONT fixture, whose base qualities run Q2-Q19. The SNP-proximity scan
    # below does still want the threshold, and applies it in column_bases().
    for pcol in bamfile.pileup(reference=chrom, start=mutstart-1, end=mutend+1, max_depth=int(args.maxdepth), ignore_overlaps=False, min_base_quality=0):
        # `if pcol.pos:` used to guard this loop, silently skipping position 0
        # of a contig.
        if args.ignorepileup and (pcol.pos < mutstart-1 or pcol.pos > mutend+1):
            continue

        is_mutcol = pcol.pos+1 in mutpos_list

        for pread in pcol.pileups:
            if avoid is not None and pread.alignment.qname in avoid:
                logger.warning(region + " dropped mutation due to read in --avoidlist " + pread.alignment.qname)
                return MutateResult(failed=True, maxfrac=result.maxfrac)

            if not is_mutcol:
                continue

            # only consider primary alignments
            if pread.query_position is not None and not pread.alignment.is_secondary and bin(pread.alignment.flag & 2048) != bin(2048):
                pairname = 'F' # read is first in pair
                if pread.alignment.is_read2:
                    pairname = 'S' # read is second in pair
                if not pread.alignment.is_paired:
                    pairname = 'U' # read is unpaired

                extqname = ','.join((pread.alignment.qname,str(pread.alignment.pos),pairname))

                if not pread.alignment.mate_is_unmapped:
                    result.outreads[extqname] = pread.alignment
                    mutid = mutid_list[mutpos_list.index(pcol.pos+1)]

                    if is_snv:
                        if extqname not in result.mutreads:
                            result.mutreads[extqname] = pread.alignment.seq

                        mutbase = mutbase_list[mutpos_list.index(pcol.pos+1)]
                        mutbases = list(result.mutreads[extqname])
                        mutbases[pread.query_position] = mutbase
                        result.mutreads[extqname] = ''.join(mutbases)

                    if is_insertion:
                        result.mutreads[extqname] = makeins(pread.alignment, indel_start, ins_seq)

                    if is_deletion:
                        result.mutreads[extqname] = makedel(pread.alignment, chrom, indel_start, indel_end, reffile)

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
                                return MutateResult(failed=True, maxfrac=result.maxfrac)

                    if extqname not in result.mutmates:
                        result.mutmates[extqname] = mate

                    log.write(" ".join(('read',extqname,result.mutreads[extqname],"\n")))

                if len(result.mutreads) > int(args.maxdepth):
                    logger.warning("depth at site is greater than cutoff, aborting mutation")
                    return MutateResult(failed=True, maxfrac=result.maxfrac)

        # make sure region doesn't have any changes that are likely SNPs
        # (trying to avoid messing with haplotypes)
        #
        # maxfrac and hasSNP used to be reset at the top of this block, so only
        # the last column examined could ever set them.
        if result.maxfrac is None:
            result.maxfrac = 0.0

        if ignoresnps:
            continue

        basepile = column_bases(pcol)
        if basepile:
            majb = majorbase(basepile)
            minb = minorbase(basepile)

            frac = float(minb[1])/(float(majb[1])+float(minb[1]))
            if minb[0] == majb[0]:
                frac = 0.0
            if frac > result.maxfrac:
                result.maxfrac = frac
            if frac > snvfrac:
                # this warning used to concatenate args.snvfrac, an int by
                # default, onto a str and raise TypeError instead of warning
                logger.warning("%s dropped for proximity to SNP, nearby SNP MAF: %f (max snv frac: %f)" % (region, frac, snvfrac))
                result.has_snp = True
        else:
            logger.warning(region + " could not pileup for region: " + chrom + ":" + str(pcol.pos))
            if not args.ignorepileup:
                result.has_snp = True

    if result.maxfrac is None:
        logger.warning("could not pile up over region: %s" % region)
        return MutateResult(failed=True)

    return result
