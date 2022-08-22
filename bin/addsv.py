#!/usr/bin/env python

#from __future__ import print_function

import re
import os
import sys
import random
import subprocess
import argparse
import pysam
import bamsurgeon.replace_reads as rr
import bamsurgeon.asmregion as ar
import bamsurgeon.mutableseq as ms
import bamsurgeon.aligners as aligners
import bamsurgeon.makevcf as makevcf

from bamsurgeon.common import *
from uuid import uuid4
from shutil import move
from collections import defaultdict as dd
from concurrent.futures import ProcessPoolExecutor

import logging
FORMAT = '%(levelname)s %(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def get_reads(bam_file, chrom, start, end, svfrac):
    bam = pysam.AlignmentFile(bam_file)
    for read in bam.fetch(chrom, start, end):
        read_end = read.reference_start + read.query_length
        pair_end = read.next_reference_start + read.query_length
        if read.is_duplicate or read.is_secondary or read.is_supplementary or \
                read.is_unmapped or read.mate_is_unmapped or \
                read.next_reference_name != chrom or \
                pair_end > end or read.next_reference_start < start or \
                read_end > end or read.reference_start < start:
            continue
        
        read_random_factor = read_hash_fraction(read.query_name)
        if read_random_factor <= svfrac:
            yield read
    bam.close()


def runwgsim(contig, newseq, pemean, pesd, tmpdir, nsimreads, mutid='null', err_rate=0.0, seed=None, trn_contig=None, rename=True):
    ''' wrapper function for wgsim, could swap out to support other reads simulators (future work?) '''

    basefn = tmpdir + '/' + mutid + ".wgsimtmp." + str(uuid4())
    fasta = basefn + ".fasta"
    fq1 = basefn + ".1.fq"
    fq2 = basefn + ".2.fq"

    fout = open(fasta,'w')
    fout.write(">target\n" + newseq + "\n")
    fout.close()

    ctg_len = len(contig)
    if trn_contig: ctg_len += len(trn_contig)

    # # adjustment factor for length of new contig vs. old contig
    logger.info("%s old ctg len: %d" % (mutid, ctg_len))
    logger.info("%s new ctg len: %d" % (mutid, len(newseq)))
    logger.info("%s num. sim. reads: %d" % (mutid, nsimreads))
    logger.info("%s PE mean outer distance: %f" % (mutid, pemean))
    logger.info("%s PE outer distance SD: %f" % (mutid, pesd))
    logger.info("%s error rate: %f" % (mutid, err_rate))

    rquals = contig.rquals
    mquals = contig.mquals

    if trn_contig:
        rquals += trn_contig.rquals
        mquals += trn_contig.mquals

    # length of quality score comes from original read, used here to set length of read
    maxqlen = 0
    for qual in (rquals + mquals):
        if len(qual) > maxqlen:
            maxqlen = len(qual)

    wgsim_args = ['wgsim','-e', str(err_rate),'-d',str(pemean),'-s',str(pesd),'-N',str(nsimreads),'-1',str(maxqlen),'-2', str(maxqlen),'-r','0','-R','0',fasta,fq1,fq2]

    seed = 1 if seed == 0 else seed # Fix for wgsim thinking 0 is no seed
    if seed is not None: wgsim_args += ['-S', str(seed)]

    logger.info(str(wgsim_args))
    subprocess.check_call(wgsim_args)

    os.remove(fasta)


    return (fq1,fq2)


def singleseqfa(file,mutid='null'):
    with open(file, 'r') as fasta:
        header = None
        seq = ''
        for line in fasta:
            line = line.strip()
            if line.startswith('>'):
                if header is not None:
                    logger.warning("%s multiple entries found in %s only using the first" % (mutid, file))
                header = line.lstrip('>')
            else:
                seq += line
    return seq


def load_inslib(infa):
    seqdict = {}

    with open(infa, 'r') as fa:
        seqid = ''
        seq   = ''
        for line in fa:
            if line.startswith('>'):
                if seq != '':
                    seqdict[seqid] = seq
                seqid = line.lstrip('>').strip()
                seq   = ''
            else:
                assert seqid != ''
                seq = seq + line.strip()

    if seqid not in seqdict and seq != '':
        seqdict[seqid] = seq

    return seqdict



def align(qryseq, refseq):
    rnd = str(uuid4())
    tgtfa = 'tmp.' + rnd + '.tgt.fa'
    qryfa = 'tmp.' + rnd + '.qry.fa'

    tgt = open(tgtfa, 'w')
    qry = open(qryfa, 'w')

    tgt.write('>ref' + '\n' + refseq + '\n')
    qry.write('>qry' + '\n' + qryseq + '\n')

    tgt.close()
    qry.close()

    cmd = ['exonerate', '--bestn', '1', '-m', 'ungapped', '--showalignment','0', '--ryo', 'SUMMARY\t%s\t%qab\t%qae\t%tab\t%tae\n', qryfa, tgtfa]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    best = []
    topscore = 0

    for pline in p.stdout.readlines():
        pline = pline.decode()
        if pline.startswith('SUMMARY'):
            c = pline.strip().split()
            if int(c[1]) > topscore:
                topscore = int(c[1])
                best = c

    os.remove(tgtfa)
    os.remove(qryfa)

    return best


def discordant_fraction(bamfile, chrom, start, end):
    r = 0
    d = 0
    bam = pysam.AlignmentFile(bamfile)
    for read in bam.fetch(chrom, start, end):
        r += 1
        if not read.is_proper_pair:
            d += 1

    if r > 0:
        return float(d)/float(r)
    else:
        return 0.0


def trim_contig(mutid, chrom, start, end, contig, reffile):
    # trim contig to get best ungapped aligned region to ref.

    refseq = reffile.fetch(chrom,start,end)
    alignstats = align(contig.seq, refseq)

    if len(alignstats) < 6:
        logger.warning("%s alignstats: %s" % (mutid, str(alignstats)))
        logger.warning("%s No good alignment between mutated contig and original, aborting mutation!" % mutid)
        return [None] * 9
    
    qrystart, qryend = map(int, alignstats[2:4])
    tgtstart, tgtend = map(int, alignstats[4:6])

    refseq = refseq[tgtstart:tgtend]
    
    logger.info("%s alignment result: %s" % (mutid, str(alignstats)))

    contig.trimseq(qrystart, qryend)
    logger.info("%s trimmed contig length: %d" % (mutid, contig.len))

    if tgtstart > tgtend: # detect reverse complemented contig
        contig.rc = True

    refstart = start + tgtstart
    refend = start + tgtend

    if refstart > refend:
        refstart, refend = refend, refstart


    return contig, refseq, alignstats, refstart, refend, qrystart, qryend, tgtstart, tgtend

def locate_contig_pos(refstart, refend, user_start, user_end, contig_len, maxlibsize):
    contig_start = None
    contig_end = None

    if user_start - refstart > maxlibsize:
        contig_start = (user_start - refstart)

    if refend - user_end > maxlibsize:
        contig_end = contig_len - (refend - user_end)

    return contig_start, contig_end



def add_donor_reads(args, mutid, tmpbamfn, bdup_chrom, bdup_left_bnd, bdup_right_bnd, bdup_svfrac):
    tmpbam = pysam.AlignmentFile(tmpbamfn)

    outbamfn = '%s/%s.%s.bigdup.merged.bam' % (args.tmpdir, mutid, str(uuid4()))
    outbam = pysam.AlignmentFile(outbamfn, 'wb', template=tmpbam)
    for read in tmpbam.fetch(until_eof=True):
        outbam.write(read)

    # Calculate donor norm factor
    with pysam.AlignmentFile(args.donorbam) as donorbam:
        cover_donor = donorbam.count(contig=bdup_chrom, start=bdup_left_bnd, end=bdup_right_bnd) / float(bdup_right_bnd-bdup_left_bnd)
    with pysam.AlignmentFile(args.bamFileName) as origbam:
        cover_orig = origbam.count(contig=bdup_chrom, start=bdup_left_bnd, end=bdup_right_bnd) / float(bdup_right_bnd-bdup_left_bnd)

    donor_norm_factor = cover_orig * bdup_svfrac / cover_donor
    if donor_norm_factor > 1.0:
        logger.warning('%s: donor_norm_factor %f > 1.0. This means donor bam has less coverage than required.' % (mutid, donor_norm_factor))

    logger.info('%s: BIGDUP donor coverage normalisation factor: %f' % (mutid, donor_norm_factor))

    logger.info('%s: fetch donor reads from %s-%d-%d' % (mutid, bdup_chrom, bdup_left_bnd, bdup_right_bnd))

    nreads = 0

    for read in get_reads(args.donorbam, bdup_chrom, bdup_left_bnd, bdup_right_bnd, donor_norm_factor):
        read.query_name = read.query_name + '_donor_' + mutid
        outbam.write(read)
        nreads += 1

    outbam.close()

    logger.info('%s: using %d donor reads from %s' % (mutid, nreads, args.donorbam))

    return outbamfn

def merge_multi_trn(args, alignopts, pair, chrom, start, end, vaf):
    assert len(pair) == 2

    mutid = os.path.basename(pair[0]).split('.')[0]

    outbamfn = '%s/%s.%s.merged.bam' % (args.tmpdir, mutid, str(uuid4()))
    bams = [pysam.AlignmentFile(bam) for bam in pair]
    outbam = pysam.AlignmentFile(outbamfn, 'wb', template=bams[0])

    readbins = {} # randomly assorted reads into bam sources 0 and 1

    for bam in bams:
        for read in bam.fetch(until_eof=True):
            readbins[read.query_name] = random.choice([0,1])

        bam.close()

    bams = [pysam.AlignmentFile(bam) for bam in pair]

    for i, bam in enumerate(bams):
        for read in bam.fetch(until_eof=True):
            if readbins[read.query_name] == i:
                outbam.write(read)

    outbam.close()

    # cleanup
    for fn in pair:
        os.remove(fn)

    return outbamfn


def makemut(args, bedline, alignopts):
    bedline = bedline.strip()

    if args.seed is not None: random.seed(args.seed + int(bedline.split()[1]))

    mutid = '_'.join(map(str, bedline.split()[:4]))

    bamfile = pysam.AlignmentFile(args.bamFileName)
    reffile = pysam.Fastafile(args.refFasta)
    logfn = '_'.join(map(os.path.basename, bedline.split()[:4])) + ".log"
    logfile = open('addsv_logs_' + os.path.basename(args.outBamFile) + '/' + os.path.basename(args.outBamFile) + '_' + logfn, 'w')
    mutinfo = {}

    # optional CNV file
    cnv = None
    if (args.cnvfile):
        cnv = pysam.Tabixfile(args.cnvfile, 'r')

    # temporary file to hold mutated reads
    outbam_mutsfile = args.tmpdir + '/' + '.'.join((mutid, str(uuid4()), "muts.bam"))

    c = bedline.split()
    chrom  = c[0]
    start  = int(c[1])
    end    = int(c[2])
    araw   = c[3:] # INV, DEL, INS,  DUP, TRN

    # desired start/end
    user_start = start
    user_end   = end

    # Check if has sufficient depth
    user_start_depth = bamfile.count(chrom, user_start-1, user_start)
    user_end_depth   = bamfile.count(chrom, user_end-1, user_end)
    if user_start_depth < args.mindepth or user_end_depth < args.mindepth:
        logger.warning('%s skipping due to insufficient depth %d %d' % (mutid, user_start_depth, user_end_depth))
        return None, None, None
    elif user_start_depth > args.maxdepth or user_end_depth > args.maxdepth:
        logger.warning('%s skipping due to excessive depth %d %d' % (mutid, user_start_depth, user_end_depth))
        return None, None, None

    # translocation specific
    trn_chrom = None
    trn_start = None
    trn_end   = None

    is_transloc = c[3] in ('TRN', 'BIGDEL', 'BIGINV', 'BIGDUP')

    if is_transloc:
        araw = [c[3]]
        if len(c) > 7:
            araw += c[7:]

        start -= int(args.minctglen)
        end   += int(args.minctglen)
        if start < 0: start = 0

        trn_chrom = c[4]
        user_trn_start = int(c[5])
        user_trn_end   = int(c[6])

        # Check for sufficient depth
        user_trn_start_depth = bamfile.count(trn_chrom, user_trn_start-1, user_trn_start)
        user_trn_end_depth   = bamfile.count(trn_chrom, user_trn_end-1, user_trn_end)
        if user_trn_start_depth < args.mindepth or user_trn_end_depth < args.mindepth:
            logger.warning('%s skipping due to insufficient depth %d %d' % (mutid, user_trn_start_depth, user_trn_end_depth))
            return None, None, None
        elif user_trn_start_depth > args.maxdepth or user_trn_end_depth > args.maxdepth:
            logger.warning('%s skipping due to excessive depth %d %d' % (mutid, user_trn_start_depth, user_trn_end_depth))
            return None, None, None

        trn_start = int(c[5]) - int(args.minctglen)
        trn_end   = int(c[6]) + int(args.minctglen)
        if trn_start < 0: trn_start = 0

    actions = map(lambda x: x.strip(),' '.join(araw).split(';'))

    svfrac = float(args.svfrac) # default, can be overridden by cnv file or per-variant

    cn = 1.0

    trn_left_flip  = False
    trn_right_flip = False

    if cnv: # CNV file is present
        if chrom in cnv.contigs:
            for cnregion in cnv.fetch(chrom,start,end):
                cn = float(cnregion.strip().split()[3]) # expect chrom,start,end,CN
                logger.info("INFO" + mutid + "\t" + ' '.join(("copy number in sv region:",chrom,str(start),str(end),"=",str(cn))) + "\n")
                svfrac = svfrac/float(cn)
                assert svfrac <= 1.0, 'copy number from %s must be at least 1: %s' % (args.cnvfile, cnregion.strip())
                logger.info("INFO" + mutid + "\tadjusted default MAF: " + str(svfrac) + "\n")

    logger.info("%s interval: %s" % (mutid, bedline))
    logger.info("%s length: %d" % (mutid, (end-start)))

   # modify start and end if interval is too short
    minctglen = int(args.minctglen)

    # adjust if minctglen is too short
    if minctglen < 3*int(args.maxlibsize):
        minctglen = 3*int(args.maxlibsize)

    #if end-start < minctglen:
    adj   = minctglen - (end-start)
    start = int(start - adj/2)
    end   = int(end + adj/2)

    #logger.info("%s note: interval size was too short, adjusted: %s:%d-%d" % (mutid, chrom, start, end))

    dfrac = discordant_fraction(args.bamFileName, chrom, start, end)
    logger.info("%s discordant fraction: %f" % (mutid, dfrac))

    if dfrac > args.maxdfrac:
        logger.warning("%s discordant fraction %f > %f aborting mutation!\n" % (mutid, dfrac, args.maxdfrac))
        return None, None, None

    contigs = ar.asm(chrom, start, end, args.bamFileName, reffile, int(args.kmersize), args.tmpdir, mutid=mutid, debug=args.debug)

    if len(contigs) == 0:
        logger.warning("%s generated no contigs, skipping site." % mutid)
        return None, None, None

    trn_contigs = None
    if is_transloc:
        logger.info("%s assemble translocation end: %s:%d-%d" % (mutid, trn_chrom, trn_start, trn_end))
        trn_contigs = ar.asm(trn_chrom, trn_start, trn_end, args.bamFileName, reffile, int(args.kmersize), args.tmpdir, mutid=mutid, debug=args.debug)

    maxcontig = sorted(contigs)[-1]

    trn_maxcontig = None
    rename_reads = True

    if is_transloc:
        if len(trn_contigs) == 0:
            logger.warning("%s translocation partner generated no contigs, skipping site." % mutid)
            return None, None, None

        trn_maxcontig = sorted(trn_contigs)[-1]

    if re.search('N', maxcontig.seq):
        if args.allowN:
            logger.warning("%s contig has ambiguous base (N), replaced with 'A'" % mutid)
            maxcontig.seq = re.sub('N', 'A', maxcontig.seq)
        else:
            logger.warning("%s contig dropped due to ambiguous base (N), aborting mutation." % mutid)
            return None, None, None

    if is_transloc and re.search('N', trn_maxcontig.seq):
        if args.allowN:
            logger.warning("%s contig has ambiguous base (N), replaced with 'A'" % mutid)
            trn_maxcontig.seq = re.sub('N', 'A', trn_maxcontig.seq)
        else:
            logger.warning("%s contig dropped due to ambiguous base (N), aborting mutation." % mutid)
            return None, None, None

    if maxcontig is None:
        logger.warning("%s maxcontig has length 0, aborting mutation!" % mutid)
        return None, None, None

    if is_transloc and trn_maxcontig is None:
        logger.warning("%s transloc maxcontig has length 0, aborting mutation!" % mutid)
        return None, None, None

    logger.info("%s best contig length: %d" % (mutid, sorted(contigs)[-1].len))

    if is_transloc:
        logger.info("%s best transloc contig length: %d" % (mutid, sorted(trn_contigs)[-1].len))

    # trim contig to get best ungapped aligned region to ref.
    maxcontig, refseq, alignstats, refstart, refend, qrystart, qryend, tgtstart, tgtend = trim_contig(mutid, chrom, start, end, maxcontig, reffile)

    if maxcontig is None:
        logger.warning("%s best contig did not have sufficent match to reference, aborting mutation." % mutid)
        return None, None, None

    logger.info("%s start: %d, end: %d, tgtstart: %d, tgtend: %d, refstart: %d, refend: %d" % (mutid, start, end, tgtstart, tgtend, refstart, refend))

    if is_transloc:
        trn_maxcontig, trn_refseq, trn_alignstats, trn_refstart, trn_refend, trn_qrystart, trn_qryend, trn_tgtstart, trn_tgtend = trim_contig(mutid, trn_chrom, trn_start, trn_end, trn_maxcontig, reffile)

        if trn_maxcontig is None:
            logger.warning("%s best contig for translocation partner did not have sufficent match to reference, aborting mutation." % mutid)
            return None, None, None
            
        logger.info("%s trn_start: %d, trn_end: %d, trn_tgtstart: %d, trn_tgtend:%d , trn_refstart: %d, trn_refend: %d" % (mutid, trn_start, trn_end, trn_tgtstart, trn_tgtend, trn_refstart, trn_refend))

    # is there anough room to make mutations?
    if maxcontig.len < 3*int(args.maxlibsize):
        logger.warning("%s best contig too short to make mutation!" % mutid)
        return None, None, None

    if is_transloc and trn_maxcontig.len < 3*int(args.maxlibsize):
        logger.warning("%s best transloc contig too short to make mutation!" % mutid)
        return None, None, None

    # make mutation in the largest contig
    mutseq = ms.MutableSeq(maxcontig.seq)

    if maxcontig.rc:
        mutseq = ms.MutableSeq(rc(maxcontig.seq)) 

    trn_mutseq = None

    if is_transloc:
        if trn_maxcontig.rc:
            trn_mutseq = ms.MutableSeq(rc(trn_maxcontig.seq))
        else:
            trn_mutseq = ms.MutableSeq(trn_maxcontig.seq)


    # support for multiple mutations
    for actionstr in actions:
        a = actionstr.split()
        action = a[0]

        logger.info("%s action: %s %s" % (mutid, actionstr, action))

        insseqfile = None
        insseq = ''
        tsdlen = 0  # target site duplication length
        ndups = 0   # number of tandem dups
        dsize = 0.0 # deletion size fraction
        dlen = 0
        ins_motif = None

        if action == 'INS':
            assert len(a) > 1 # insertion syntax: INS <file.fa> [optional TSDlen]
            insseqfile = a[1]
            if not (os.path.exists(insseqfile) or insseqfile == 'RND' or insseqfile.startswith('INSLIB:')): # not a file... is it a sequence? (support indel ins.)
                assert re.search('^[ATGCatgc]*$',insseqfile), "cannot determine SV type: %s" % insseqfile # make sure it's a sequence
                insseq = insseqfile.upper()
                insseqfile = None
            if len(a) > 2: # field 5 for insertion is TSD Length
                tsdlen = int(a[2])

            if len(a) > 3: # field 6 for insertion is motif, format = 'NNNN^NNNN where ^ is cut site
                ins_motif = a[3]
                assert '^' in ins_motif, 'insertion motif specification requires cut site defined by ^'

            if len(a) > 4: # field 7 is VAF
                svfrac = float(a[4])/cn

        if action == 'DUP':
            if len(a) > 1:
                ndups = int(a[1])
            else:
                ndups = 1

            if len(a) > 2: # VAF
                svfrac = float(a[2])/cn

        if action == 'DEL':
            dsize = 1.0

            if len(a) > 1: # VAF
                svfrac = float(a[1])/cn

        if action in ('TRN', 'BIGDEL', 'BIGINV', 'BIGDUP'):
            if len(a) > 1: # translocation end orientation ++ / +- / -+ / --
                trn_left_flip = a[1][0] == '-'
                trn_right_flip = a[1][1] == '-'

            if len(a) > 2:
                svfrac = float(a[2])/cn

        if action == 'INV':
            if len(a) > 1:
                svfrac = float(a[1])/cn


        logger.info("%s final VAF accounting for copy number %f: %f" % (mutid, cn, svfrac))

        logfile.write(">" + chrom + ":" + str(refstart) + "-" + str(refend) + " BEFORE\n" + str(mutseq) + "\n")

        contig_start = None
        contig_end = None
        trn_contig_start = None
        trn_contig_end = None
        exact_success = True

        contig_start, contig_end = locate_contig_pos(refstart, refend, user_start, user_end, mutseq.length(), int(args.maxlibsize))

        if contig_start is None:
            logger.warning('%s contig does not cover user start' % mutid)
            exact_success = False
            #print refstart, refend, user_start, user_end, int(args.maxlibsize)

        if contig_end is None:
            logger.warning('%s contig does not cover user end' % mutid)
            exact_success = False
            #print refstart, refend, user_start, user_end, int(args.maxlibsize)

        if is_transloc:
            trn_contig_start, trn_contig_end = locate_contig_pos(trn_refstart, trn_refend, user_trn_start, user_trn_end, trn_mutseq.length(), int(args.maxlibsize))

            if trn_contig_start is None:
                logger.warning('%s contig does not cover user translocation start' % mutid)
                exact_success = False

            if trn_contig_end is None:
                logger.warning('%s contig does not cover user translocation end' % mutid)
                exact_success = False


        if args.require_exact and not exact_success:
            logger.warning('%s dropped mutation due to --require_exact')
            return None, None, None


        if action == 'INS':
            inspoint = int(mutseq.length()/2)
            if None not in (contig_start, contig_end):
                inspoint = int((contig_start+contig_end)/2)

            if ins_motif is not None:
                inspoint = mutseq.find_site(ins_motif, left_trim=int(args.maxlibsize), right_trim=int(args.maxlibsize))

                if inspoint < int(args.maxlibsize) or inspoint > mutseq.length() - int(args.maxlibsize):
                    logger.info("%s picked midpoint, no cutsite found" % mutid)
                    inspoint = int(mutseq.length()/2)

            if insseqfile: # seq in file
                if insseqfile == 'RND':
                    assert args.inslib is not None # insertion library needs to exist
                    insseqfile = random.choice(list(args.inslib.keys()))
                    logger.info("%s chose sequence from insertion library: %s" % (mutid, insseqfile))
                    mutseq.insertion(inspoint, args.inslib[insseqfile], tsdlen)

                elif insseqfile.startswith('INSLIB:'):
                    assert args.inslib is not None # insertion library needs to exist
                    insseqfile = insseqfile.split(':')[1]
                    logger.info("%s specify sequence from insertion library: %s" % (mutid, insseqfile))
                    assert insseqfile in args.inslib, '%s not found in insertion library' % insseqfile
                    mutseq.insertion(inspoint, args.inslib[insseqfile], tsdlen)

                else:
                    mutseq.insertion(inspoint, singleseqfa(insseqfile, mutid=mutid), tsdlen)

            else: # seq is input
                mutseq.insertion(inspoint, insseq, tsdlen)

            ins_len = len(mutseq.seq) - len(maxcontig.seq)
            mutinfo[mutid] = "\t".join(('ins',chrom,str(refstart),str(refend),action,str(mutseq.length()),str(inspoint),str(insseqfile),str(tsdlen),str(ins_len),str(svfrac)))
            logfile.write(mutinfo[mutid] + "\n")

        elif action == 'INV':
            invstart = int(args.maxlibsize)
            invend = mutseq.length() - invstart

            if None not in (contig_start, contig_end):
                invstart = contig_start
                invend   = contig_end

            mutseq.inversion(invstart,invend)

            mutinfo[mutid] = "\t".join(('inv',chrom,str(refstart),str(refend),action,str(mutseq.length()),str(invstart),str(invend),str(svfrac)))
            logfile.write(mutinfo[mutid] + "\n")

        elif action == 'DEL':
            delstart = int(args.maxlibsize)
            delend = mutseq.length() - delstart

            if None not in (contig_start, contig_end):
                delstart = contig_start
                delend   = contig_end

            mutseq.deletion(delstart,delend)

            mutinfo[mutid] = "\t".join(('del',chrom,str(refstart),str(refend),action,str(mutseq.length()),str(delstart),str(delend),str(dlen),str(svfrac)))
            logfile.write(mutinfo[mutid] + "\n")

        elif action == 'DUP':
            dupstart = int(args.maxlibsize)
            dupend = mutseq.length() - dupstart

            if None not in (contig_start, contig_end):
                dupstart = contig_start
                dupend   = contig_end

            mutseq.duplication(dupstart,dupend,ndups)

            mutinfo[mutid] = "\t".join(('dup',chrom,str(refstart),str(refend),action,str(mutseq.length()),str(dupstart),str(dupend),str(ndups),str(svfrac)))
            logfile.write(mutinfo[mutid] + "\n")

        elif action == 'TRN':
            trnpoint_1 = int(mutseq.length()/2)
            trnpoint_2 = int(trn_mutseq.length()/2)

            if None not in (contig_start, contig_end):
                trnpoint_1 = int((contig_start + contig_end)/2)

            if None not in (trn_contig_start, trn_contig_end):
                trnpoint_2 = int((trn_contig_start + trn_contig_end)/2)

            mutseq.fusion(trnpoint_1, trn_mutseq, trnpoint_2, flip1=trn_left_flip, flip2=trn_right_flip)

            mutinfo[mutid] = "\t".join(('trn',chrom,str(refstart),str(refend),action,str(trnpoint_1),trn_chrom,str(trn_refstart),str(trn_refend),str(trnpoint_2),str(trn_left_flip),str(trn_right_flip),str(svfrac)))
            logfile.write(mutinfo[mutid] + "\n")

        elif action == 'BIGDEL':
            trnpoint_1 = int(mutseq.length()/2)
            trnpoint_2 = int(trn_mutseq.length()/2)

            if None not in (contig_start, contig_end):
                trnpoint_1 = int((contig_start + contig_end)/2)

            if None not in (trn_contig_start, trn_contig_end):
                trnpoint_2 = int((trn_contig_start + trn_contig_end)/2)

            mutseq.fusion(trnpoint_1, trn_mutseq, trnpoint_2)

            mutinfo[mutid] = "\t".join(('bigdel',chrom,str(refstart),str(refend),action,str(trnpoint_1),trn_chrom,str(trn_refstart),str(trn_refend),str(trnpoint_2),str(svfrac)))
            logfile.write(mutinfo[mutid] + "\n")

        elif action == 'BIGINV':
            trnpoint_1 = int(mutseq.length()/2)
            trnpoint_2 = int(trn_mutseq.length()/2)

            if None not in (contig_start, contig_end):
                trnpoint_1 = int((contig_start + contig_end)/2)

            if None not in (trn_contig_start, trn_contig_end):
                trnpoint_2 = int((trn_contig_start + trn_contig_end)/2)

            mutseq.fusion(trnpoint_1, trn_mutseq, trnpoint_2, flip1=trn_left_flip, flip2=trn_right_flip)

            mutinfo[mutid] = "\t".join(('biginv',chrom,str(refstart),str(refend),action,str(trnpoint_1),trn_chrom,str(trn_refstart),str(trn_refend),str(trnpoint_2),str(svfrac)))
            logfile.write(mutinfo[mutid] + "\n")

        elif action == 'BIGDUP':
            trnpoint_1 = int(mutseq.length()/2)
            trnpoint_2 = int(trn_mutseq.length()/2)

            if None not in (contig_start, contig_end):
                trnpoint_1 = int((contig_start + contig_end)/2)

            if None not in (trn_contig_start, trn_contig_end):
                trnpoint_2 = int((trn_contig_start + trn_contig_end)/2)

            mutseq.fusion(trnpoint_1, trn_mutseq, trnpoint_2)

            mutinfo[mutid] = "\t".join(('bigdup',chrom,str(refstart),str(refend),action,str(trnpoint_1),trn_chrom,str(trn_refstart),str(trn_refend),str(trnpoint_2),str(svfrac)))
            logfile.write(mutinfo[mutid] + "\n")
            rename_reads = False

        else:
            raise ValueError("ERROR " + mutid + "\t: mutation not one of: INS,INV,DEL,DUP,TRN,BIGDEL,BIGINV,BIGDUP\n")

        logfile.write(">" + chrom + ":" + str(refstart) + "-" + str(refend) +" AFTER\n" + str(mutseq) + "\n")

    pemean, pesd = float(args.ismean), float(args.issd) 
    logger.info("%s set paired end mean distance: %f" % (mutid, pemean))
    logger.info("%s set paired end distance stddev: %f" % (mutid, pesd))

    exclfile = args.tmpdir + '/' + '.'.join((mutid, 'exclude', str(uuid4()), 'txt'))
    exclude = open(exclfile, 'w')

    if is_transloc:
        buffer = int(float(args.ismean))
        region_1_start, region_1_end = (refstart + trnpoint_1 - buffer, refend) if trn_left_flip else (refstart, refstart + trnpoint_1 + buffer)
        region_2_start, region_2_end = (trn_refstart + trnpoint_2 - buffer, trn_refend) if not trn_right_flip else (trn_refstart, trn_refstart + trnpoint_2 + buffer)
        region_1_reads = get_reads(args.bamFileName, chrom, region_1_start, region_1_end, float(svfrac))
        region_2_reads = get_reads(args.bamFileName, trn_chrom, region_2_start, region_2_end, float(svfrac))
        excl_reads_names = set([read.query_name for read in region_1_reads] + [read.query_name for read in region_2_reads]) 
        nsimreads = len(excl_reads_names)
        # add additional excluded reads if bigdel(s) present
        if action == 'BIGDEL':
            bigdel_region_reads = get_reads(args.bamFileName, chrom, region_1_start, region_2_end, float(svfrac))
            excl_reads_names = set([read.query_name for read in bigdel_region_reads])
    else:
        region_reads = get_reads(args.bamFileName, chrom, refstart, refend, float(svfrac))
        excl_reads_names = set([read.query_name for read in region_reads])
        reads_ratio = len(mutseq.seq) / len(maxcontig.seq)
        nsimreads = int(len(excl_reads_names) * reads_ratio)

    for name in excl_reads_names:
        exclude.write(name + "\n")
    exclude.close()

    # simulate reads
    (fq1, fq2) = runwgsim(maxcontig, mutseq.seq, pemean, pesd, args.tmpdir, nsimreads, err_rate=float(args.simerr), mutid=mutid, seed=args.seed, trn_contig=trn_maxcontig, rename=rename_reads)

    outreads = aligners.remap_fastq(args.aligner, fq1, fq2, args.refFasta, outbam_mutsfile, alignopts, mutid=mutid, threads=int(args.alignerthreads))

    if outreads == 0:
        logger.warning("%s outbam %s has no mapped reads!" % (mutid, outbam_mutsfile))
        # Remove content from logfile in order to skip this mutation in the final VCF file
        logfile.seek(0)
        logfile.truncate()
        return None, None, None

    if action == 'BIGDUP':
        bdup_left_bnd = min(region_1_start, region_2_start, region_1_end, region_2_end)
        bdup_right_bnd = max(region_1_start, region_2_start, region_1_end, region_2_end)
        prev_outbam_mutsfile = outbam_mutsfile
        outbam_mutsfile = add_donor_reads(args, mutid, outbam_mutsfile, chrom, bdup_left_bnd, bdup_right_bnd, float(svfrac))
        os.remove(prev_outbam_mutsfile)
        os.remove(prev_outbam_mutsfile + '.bai')

    logger.info("%s temporary bam: %s" % (mutid, outbam_mutsfile))

    bamfile.close()

    return outbam_mutsfile, exclfile, mutinfo




def main(args):
    logger.info("starting %s called with args: %s" % (sys.argv[0], ' '.join(sys.argv)))
    tmpbams = [] # temporary BAMs, each holds the realigned reads for one mutation
    exclfns = [] # 'exclude' files store reads to be removed from the original BAM due to deletions

    if (args.bamFileName.endswith('.bam') and not os.path.exists(args.bamFileName + '.bai')) or \
        (args.bamFileName.endswith('.cram') and not os.path.exists(args.bamFileName + '.crai')):
        logger.error("input file must be indexed, not .bai or .crai file found for %s" % args.bamFileName)
        sys.exit(1)

    alignopts = {}
    if args.alignopts is not None:
        alignopts = dict([o.split(':') for o in args.alignopts.split(',')])

    aligners.checkoptions(args.aligner, alignopts, None, sv=True)

    # load insertion library if present
    try:
        if args.inslib is not None:
            logger.info("loading insertion library from %s" % args.inslib)
            args.inslib = load_inslib(args.inslib)
    except Exception:
        logger.error("failed to load insertion library %s" % args.inslib)
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)

    results = []
    pool = ProcessPoolExecutor(max_workers=int(args.procs))

    nmuts = 0

    if not os.path.exists(args.tmpdir):
        os.mkdir(args.tmpdir)
        logger.info("created tmp directory: %s" % args.tmpdir)

    if not os.path.exists('addsv_logs_' + os.path.basename(args.outBamFile)):
        os.mkdir('addsv_logs_' + os.path.basename(args.outBamFile))
        logger.info("created log directory: addsv_logs_%s" % os.path.basename(args.outBamFile))

    assert os.path.exists('addsv_logs_' + os.path.basename(args.outBamFile)), "could not create output directory!"
    assert os.path.exists(args.tmpdir), "could not create temporary directory!"

    biginvs = {}

    with open(args.varFileName, 'r') as varfile:
        for bedline in varfile:
            bedline = bedline.strip()
            multi_part = []

            if re.search('^#',bedline):
                continue

            if args.maxmuts and nmuts >= int(args.maxmuts):
                break

            mut_type = bedline.split()[3]
            mut_len = int(bedline.split()[2]) - int(bedline.split()[1])

            if mut_type in ('DEL', 'DUP', 'INV') and mut_len > 10000:
                logger.warning('%s is over 10kbp long: converting to use BIG%s instead.' % (bedline, mut_type))
                mut_type = 'BIG' + mut_type
                if mut_type == 'BIGDUP' and len(bedline.split()) == 6: # convert DUP to BIGDUP
                    b = bedline.split()
                    bedline = ' '.join((b[:4] + [b[-1]]))

            if mut_type.startswith('BIG') and mut_len < 5000:
                mut_type = mut_type.replace('BIG', '')
                logger.warning('%s is under 5kbp, "BIG" mutation types will yield unpredictable results, converting to %s' % (bedline, mut_type))

            # rewrite bigdel coords as translocation
            if mut_type == 'BIGDEL':
                bdel_svfrac = float(args.svfrac)
                if len(bedline.split()) == 5:
                    bdel_svfrac = float(bedline.split()[-1])

                bdel_chrom, bdel_start, bdel_end = bedline.split()[:3]
                bdel_start = int(bdel_start)
                bdel_end   = int(bdel_end)

                bdel_left_start = bdel_start
                bdel_left_end   = bdel_start

                bdel_right_start = bdel_end
                bdel_right_end   = bdel_end

                bedline = '%s %d %d BIGDEL %s %d %d %s %f' % (bdel_chrom, bdel_left_start, bdel_left_end, bdel_chrom, bdel_right_start, bdel_right_end, '++', bdel_svfrac)

            # rewrite bigdup coords as translocation
            if mut_type == 'BIGDUP':
                bdup_svfrac = float(args.svfrac)
                if len(bedline.split()) == 6:
                    bdup_svfrac = float(bedline.split()[-1])

                if args.donorbam is None:
                    logger.warning('%s: using BIGDUP requires specifying a --donorbam and none was provided, using %s' % (bedline, args.bamFileName))
                    args.donorbam = args.bamFileName

                bdup_chrom, bdup_start, bdup_end = bedline.split()[:3]
                bdup_start = int(bdup_start)
                bdup_end   = int(bdup_end)

                bdup_left_start = bdup_start
                bdup_left_end   = bdup_start

                bdup_right_start = bdup_end
                bdup_right_end   = bdup_end

                bedline = '%s %d %d BIGDUP %s %d %d %s %f' % (bdup_chrom, bdup_right_start, bdup_right_end, bdup_chrom, bdup_left_start, bdup_left_end, '++', bdup_svfrac)

            # rewrite biginv coords as translocations
            if mut_type == 'BIGINV':

                binv_svfrac = float(args.svfrac)
                if len(bedline.split()) == 5:
                    binv_svfrac = float(bedline.split()[-1])

                binv_chrom, binv_start, binv_end = bedline.split()[:3]
                binv_start = int(binv_start)
                binv_end = int(binv_end)

                binv_left_start = binv_start
                binv_left_end   = binv_start

                binv_right_start = binv_end
                binv_right_end   = binv_end

                # left breakpoint
                multi_part.append('%s %d %d BIGINV %s %d %d %s %f' % (binv_chrom, binv_left_start, binv_left_end, binv_chrom, binv_right_start, binv_right_end, '+-', binv_svfrac))

                # right breakpoint
                multi_part.append('%s %d %d BIGINV %s %d %d %s %f' % (binv_chrom, binv_left_start, binv_left_end, binv_chrom, binv_right_start, binv_right_end, '-+', binv_svfrac))

                binv_mutid = '_'.join(map(str, (binv_chrom, binv_left_start, binv_left_end, 'BIGINV')))
                biginvs[binv_mutid] = (binv_chrom, binv_start, binv_end, binv_svfrac)

            if len(multi_part) == 0:
                # submit each mutation as its own thread
                result = pool.submit(makemut, args, bedline, alignopts)
                results.append(result)

            else:
                for bedline in multi_part:
                    result = pool.submit(makemut, args, bedline, alignopts)
                    results.append(result)

            nmuts += 1

    ## process the results of mutation jobs
    for result in results:
        tmpbam = None
        exclfn = None

        tmpbam, exclfn, mutinfo = result.result()

        if None not in (tmpbam, exclfn) and os.path.exists(tmpbam) and os.path.exists(exclfn):
            if bamreadcount(tmpbam) > 0:
                tmpbams.append(tmpbam)
                exclfns.append(exclfn)
            else:
                os.remove(tmpbam)
                os.remove(exclfn)

    if len(tmpbams) == 0:
        logger.error("no succesful mutations")
        sys.exit(1)

    biginv_pairs = dd(list)

    new_tmpbams = []

    for tmpbamfn in tmpbams:
        mutid = os.path.basename(tmpbamfn).split('.')[0]
        if mutid.endswith('BIGINV'):
            biginv_pairs[mutid].append(tmpbamfn)
        else:
            new_tmpbams.append(tmpbamfn)

    # find translocation pairs corresponding to BIGINV, merge pairs / remove singletons
    for binv_pair in biginv_pairs.values():
        if len(binv_pair) == 2:
            logger.info('merging biginv pair and reversing unassembled interval: %s' % str(binv_pair))

            binv_mutid = os.path.basename(binv_pair[0]).split('.')[0]

            assert binv_mutid in biginvs

            binv_chrom, binv_start, binv_end, binv_svfrac = biginvs[binv_mutid]

            binv_left_end  = int(binv_left_end)
            binv_right_end = int(binv_right_end)

            if binv_left_end > binv_right_end:
                binv_left_end, binv_right_end = binv_right_end, binv_left_end

            merged_binv = merge_multi_trn(args, alignopts, binv_pair, binv_chrom, binv_start, binv_end, binv_svfrac)
            new_tmpbams.append(merged_binv)

    tmpbams = new_tmpbams

    logger.info("tmpbams: %s" % tmpbams)
    logger.info("exclude: %s" % exclfns)

    if len(tmpbams) == 0:
        sys.exit('no tmp bams remain, nothing to do!')

    excl_merged = 'addsv.exclude.final.' + str(uuid4()) + '.txt'
    mergedtmp = 'addsv.mergetmp.final.' + str(uuid4()) + '.bam'

    logger.info("merging exclude files into %s" % excl_merged)
    exclout = open(excl_merged, 'w')
    for exclfn in exclfns:
        with open(exclfn, 'r') as excl:
            for line in excl:
                exclout.write(line)
        if not args.debug:
            os.remove(exclfn)
    exclout.close()

    if len(tmpbams) == 1:
        logger.info("only one bam: %s renaming to %s" % (tmpbams[0], mergedtmp))
        os.rename(tmpbams[0], mergedtmp)
    elif len(tmpbams) > 1:
        logger.info("merging bams into %s" % mergedtmp)
        mergebams(tmpbams, mergedtmp, debug=args.debug)

    if args.skipmerge:
        logger.info("final merge skipped, please merge manually: %s" % mergedtmp)
        logger.info("exclude file to use: %s" % excl_merged)
    else:
        if args.tagreads:
            from bamsurgeon.markreads import markreads
            tmp_tag_bam = 'tag.%s.bam' % str(uuid4())
            markreads(mergedtmp, tmp_tag_bam)
            move(tmp_tag_bam, mergedtmp)
            logger.info("tagged reads.")

        logger.info("writing to %s" % args.outBamFile)
        rr.replace_reads(args.bamFileName, mergedtmp, args.outBamFile, excludefile=excl_merged, allreads=True, keepsecondary=args.keepsecondary, seed=args.seed, quiet=True)
        if not args.debug:
            os.remove(excl_merged)
            os.remove(mergedtmp)

        logger.info("done.")

    if not args.debug:
        for tmpbam in tmpbams:
            if os.path.isfile(tmpbam):
                os.remove(tmpbam)
            if os.path.isfile(tmpbam + '.bai'):
                os.remove(tmpbam + '.bai')

    var_basename = '.'.join(os.path.basename(args.varFileName).split('.')[:-1])
    bam_basename = '.'.join(os.path.basename(args.outBamFile).split('.')[:-1])

    vcf_fn = bam_basename + '.addsv.' + var_basename + '.vcf'

    makevcf.write_vcf_sv('addsv_logs_' + os.path.basename(args.outBamFile), args.refFasta, vcf_fn)

    logger.info('vcf output written to ' + vcf_fn)

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='adds SVs to reads, outputs modified reads as .bam along with mates')
    parser.add_argument('-v', '--varfile', dest='varFileName', required=True,
                        help='whitespace-delimited target regions for SV spike-in, see manual for syntax')
    parser.add_argument('-f', '--bamfile', dest='bamFileName', required=True,
                        help='sam/bam file from which to obtain reads')
    parser.add_argument('-r', '--reference', dest='refFasta', required=True,
                        help='reference genome, fasta indexed with bwa index _and_ samtools faidx')
    parser.add_argument('-o', '--outbam', dest='outBamFile', required=True,
                        help='.bam file name for output')
    parser.add_argument('-l', '--maxlibsize', dest='maxlibsize', default=600,
                        help="maximum fragment length of seq. library")
    parser.add_argument('-k', '--kmer', dest='kmersize', default=31, 
                        help="kmer size for assembly (default = 31)")
    parser.add_argument('-s', '--svfrac', dest='svfrac', default=1.0, 
                        help="allele fraction of variant (default = 1.0)")
    parser.add_argument('--require_exact', default=False, action='store_true',
                        help="drop mutation if breakpoints cannot be made exactly as input")
    parser.add_argument('--mindepth', default=10, type=int,
                        help='minimum read depth in the breakend position to make mutation (default = 10)')
    parser.add_argument('--maxdepth', default=2000, type=int,
                        help='maximum read depth in the breakend position to make mutation (default = 2000)')
    parser.add_argument('--maxdfrac', default=0.1, type=float,
                        help='maximum discordant fraction (is_proper_pair / is_pair) of reads (default = 0.1)')
    parser.add_argument('--minctglen', dest='minctglen', default=4000,
                        help="minimum length for contig generation, also used to pad assembly (default=4000)")
    parser.add_argument('-n', dest='maxmuts', default=None,
                        help="maximum number of mutations to make")
    parser.add_argument('-c', '--cnvfile', dest='cnvfile', default=None, 
                        help="tabix-indexed list of genome-wide absolute copy number values (e.g. 2 alleles = no change)")
    parser.add_argument('--donorbam', dest='donorbam', default=None,
                        help='bam file for donor reads if using BIGDUP mutations')
    parser.add_argument('--ismean', dest='ismean', default=300, 
                        help="mean insert size (default = estimate from region)")
    parser.add_argument('--issd', dest='issd', default=70, 
                        help="insert size standard deviation (default = estimate from region)")
    parser.add_argument('--simerr', dest='simerr', default=0.0,
                        help='error rate for wgsim-generated reads')
    parser.add_argument('-p', '--procs', dest='procs', default=1, 
                        help="split into multiple processes (default=1)")
    parser.add_argument('--inslib', default=None,
                        help='FASTA file containing library of possible insertions, use INS RND instead of INS filename to pick one')
    parser.add_argument('--aligner', default='backtrack',
                        help='supported aligners: ' + ','.join(aligners.supported_aligners_fastq))
    parser.add_argument('--alignopts', default=None,
                        help='aligner-specific options as comma delimited list of option1:value1,option2:value2,...')
    parser.add_argument('--alignerthreads', default=1,
                        help='threads used per realignment (default = 1)')
    parser.add_argument('--tagreads', action='store_true', default=False,
                        help='add BS tag to altered reads')
    parser.add_argument('--skipmerge', action='store_true', default=False,
                        help='do not merge spike-in reads back into original BAM')
    parser.add_argument('--keepsecondary', action='store_true', default=False,
                        help='keep secondary reads in final BAM')
    parser.add_argument('--debug', action='store_true', default=False,
                        help='output read tracking info to debug file, retain all intermediates')
    parser.add_argument('--tmpdir', default='addsv.tmp',
                        help='temporary directory (default=addsv.tmp)')
    parser.add_argument('--seed', default=None, type=int,
                        help='seed random number generation')
    parser.add_argument('--allowN', action='store_true', default=False,
                        help='allow N in contigs, replace with A and warn user (default: drop mutation)')
    args = parser.parse_args()
    main(args)

