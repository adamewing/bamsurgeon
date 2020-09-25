#!/usr/bin/env python

#from __future__ import print_function

import re
import os
import sys
import random
import subprocess
import argparse
import pysam
import bamsurgeon.replacereads as rr
import bamsurgeon.asmregion as ar
import bamsurgeon.mutableseq as ms
import bamsurgeon.aligners as aligners
import bamsurgeon.makevcf as makevcf

from bamsurgeon.common import *
from uuid import uuid4
from time import sleep
from shutil import move
from math import sqrt
from collections import Counter
from collections import defaultdict as dd
from multiprocessing import Pool

import logging
FORMAT = '%(levelname)s %(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def runwgsim(contig, newseq, svfrac, svtype, exclude, pemean, pesd, tmpdir, mutid='null', err_rate=0.0, seed=None, trn_contig=None, rename=True):
    ''' wrapper function for wgsim, could swap out to support other reads simulators (future work?) '''

    readnames = [read.name for read in contig.reads.reads.values()]
    if trn_contig: readnames += [read.name for read in trn_contig.reads.reads.values()]

    namecount = Counter(readnames)

    basefn = tmpdir + '/' + mutid + ".wgsimtmp." + str(uuid4())
    fasta = basefn + ".fasta"
    fq1 = basefn + ".1.fq"
    fq2 = basefn + ".2.fq"

    fout = open(fasta,'w')
    fout.write(">target\n" + newseq + "\n")
    fout.close()

    totalreads = len(readnames)
    paired = 0
    single = 0
    discard = 0
    pairednames = []
    # names with count 2 had both pairs in the contig
    for name,count in namecount.items():
        #print name,count
        if count == 1:
            single += 1
        elif count == 2:
            paired += 1 
            pairednames.append(name)
        else:
            discard += 1

    ctg_len = len(contig)
    if trn_contig: ctg_len += len(trn_contig)

    logger.info("%s paired reads: %d" %  (mutid, paired))
    logger.info("%s single reads: %d" % (mutid, single))
    logger.info("%s discard reads: %d" % (mutid, discard))
    logger.info("%s total reads: %d" % (mutid, totalreads))

    # adjustment factor for length of new contig vs. old contig
    lenfrac = float(len(newseq))/float(ctg_len)

    logger.info("%s old ctg len: %d" % (mutid, ctg_len))
    logger.info("%s new ctg len: %d" % (mutid, len(newseq)))
    logger.info("%s adj. factor: %f" % (mutid, lenfrac))

    # number of paried reads to simulate
    nsimreads = int((paired + (single/2)) * svfrac * lenfrac)

    logger.info("%s num. sim. reads: %d" % (mutid, nsimreads))
    logger.info("%s PE mean outer distance: %f" % (mutid, pemean))
    logger.info("%s PE outer distance SD: %f" % (mutid, pesd))
    logger.info("%s rerror rate: %f" % (mutid, err_rate))

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

    args = ['wgsim','-e', str(err_rate),'-d',str(pemean),'-s',str(pesd),'-N',str(nsimreads),'-1',str(maxqlen),'-2', str(maxqlen),'-r','0','-R','0',fasta,fq1,fq2]

    if seed is not None: args += ['-S', str(seed)]

    logger.info(str(args))
    subprocess.call(args)

    os.remove(fasta)

    if rename:
        fqReplaceList(fq1, pairednames, rquals, svfrac, svtype, exclude, mutid=mutid)
        fqReplaceList(fq2, pairednames, mquals, svfrac, svtype, exclude, mutid=mutid)

    return (fq1,fq2)


def fqReplaceList(fqfile, names, quals, svfrac, svtype, exclude, mutid='null'):
    '''
    Replace seq names in paired fastq files from a list until the list runs out
    (then stick with original names). fqfile = fastq file, names = list

    'exclude' is a filehandle, the exclude file contains read names that should
    not appear in the final output BAM

    '''
    fqin = open(fqfile,'r')

    ln = 0
    namenum = 0
    newnames = []
    seqs = []
    usednames = {}

    for fqline in fqin:
        if ln == 0:
            if len(names) > namenum:
                newnames.append(names[namenum])
            else:
                simname = fqline.strip().lstrip('@')
                simname = re.sub('/1$','',simname)  #wgsim
                simname = re.sub('/2$','',simname)  #wgsim
                newnames.append(simname) 
            namenum += 1
            ln += 1
        elif ln == 1:
            seqs.append(fqline.strip())
            ln += 1
        elif ln == 2:
            ln += 1
        elif ln == 3:
            ln = 0
        else:
            raise ValueError("ERROR\t" + now() + "\t" + mutid + "\tfastq iteration problem\n")

    fqin.close()
    os.remove(fqfile)

    # make sure there's enough (bogus) quality scores
    while len(seqs) > len(quals):
        i = random.randint(0,len(quals)-1)
        quals.append(quals[i])

    # write .fq with new names
    fqout = open(fqfile,'w')
    for i in range(namenum):
        fqout.write("@" + newnames[i] + "\n")

        # make sure quality strings are the same length as the sequences
        while len(seqs[i]) > len(quals[i]):
            quals[i] = quals[i] + 'B'

        if len(seqs[i]) < len(quals[i]):
            quals[i] = quals[i][:len(seqs[i])]

        fqout.write(seqs[i] + "\n+\n" + quals[i] + "\n")
        if newnames[i] in usednames:
            logger.warning("%s warning, used read name: %s in multiple pairs" % (mutid, newnames[i]))
        usednames[newnames[i]] = True

    is_del = False
    for sv in svtype:
        if re.search('DEL', sv):
            is_del = True

    # burn off excess if deletion
    if is_del:
        if len(seqs) > 0:
            for name in names:
                if name not in usednames:
                    if random.uniform(0,1) < svfrac:  # this controls deletion depth
                        exclude.write(name + "\n")

    fqout.close()


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


def replace(origbamfile, mutbamfile, outbamfile, excludefile, keepsecondary=False, seed=None):
    ''' open .bam file and call replacereads
    '''
    origbam = pysam.Samfile(origbamfile, 'rb')
    mutbam  = pysam.Samfile(mutbamfile, 'rb')
    outbam  = pysam.Samfile(outbamfile, 'wb', template=origbam)

    rr.replaceReads(origbam, mutbam, outbam, excludefile=excludefile, allreads=True, keepsecondary=keepsecondary, seed=seed)

    origbam.close()
    mutbam.close()
    outbam.close()


def discordant_fraction(bamfile, chrom, start, end):
    r = 0
    d = 0
    bam = pysam.Samfile(bamfile, 'rb')
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



def add_donor_reads(args, mutid, tmpbamfn, bdup_chrom, bdup_left_bnd, bdup_right_bnd, buf=200):
    assert bdup_left_bnd < bdup_right_bnd, '%s: bdup_left_bnd > bdup_right_bnd' % mutid

    donorbam = pysam.AlignmentFile(args.donorbam)

    tmpbam = pysam.AlignmentFile(tmpbamfn)

    outbamfn = '%s/%s.%s.bigdup.merged.bam' % (args.tmpdir, mutid, str(uuid4()))
    outbam = pysam.AlignmentFile(outbamfn, 'wb', template=tmpbam)

    # identify zero coverage region

    left_zero  = None
    right_zero = None

    left_cover  = None
    right_cover = None

    region = '%s:%d-%d' % (bdup_chrom, bdup_left_bnd+buf, bdup_right_bnd-buf)

    args = ['samtools', 'mpileup', '-r', region, '-a', tmpbamfn]

    FNULL = open(os.devnull, 'w')

    p = subprocess.Popen(args,stdout=subprocess.PIPE,stderr=FNULL)

    for line in p.stdout:
        line = line.decode()
        c = line.strip().split()
        pos   = int(c[1])
        depth = int(c[3])

        if left_zero is None and depth == 0:
            left_zero = pos

        if left_cover is None and depth > 0:
            left_cover = pos

        if left_zero is not None and depth == 0:
            right_zero = pos

        if depth > 0 and right_zero is not None:
            right_cover = pos

    if right_cover is None:
        right_cover = right_zero+1

    logger.info('%s: left_zero=%d, left_cover=%d, right_zero=%d, right_cover=%d' % (mutid, left_zero, left_cover, right_zero, right_cover))

    if left_cover > left_zero:
       logger.warning('%s: left_cover > left_zero' % mutid)
       left_cover, left_zero = left_zero, left_cover

    if right_cover < right_zero:
       logger.warning('%s: right_cover < right_zero' % mutid)
       right_cover, right_zero = right_zero, right_cover

    assert left_zero < right_zero, 'left_zero: %d, right_zero: %d' % (left_zero, right_zero)

    count_left  = tmpbam.count(reference=bdup_chrom, start=left_cover, end=left_zero)
    count_right = tmpbam.count(reference=bdup_chrom, start=right_zero, end=right_cover)

    cover_donor = donorbam.count(region=region) / float(bdup_right_bnd-bdup_left_bnd)

    tmpbam.reset()
    donorbam.reset()

    cover_tmp = float(count_left+count_right) / float((left_zero-left_cover)+(right_cover-right_zero))

    for read in tmpbam.fetch(until_eof=True):
        outbam.write(read)

    #donor_norm_factor = min(cover_tmp,cover_donor)/max(cover_tmp,cover_donor)
    donor_norm_factor = 1.0 # FIXME 

    logger.info('%s: BIGDUP donor coverage normalisation factor: %f' % (mutid, donor_norm_factor))

    matepairs = {}

    logger.info('%s: fetch donor reads from %s-%d-%d' % (mutid, bdup_chrom, bdup_left_bnd, bdup_right_bnd))

    paircount = 0

    for read in donorbam.fetch(bdup_chrom, bdup_left_bnd, bdup_right_bnd):
        if not read.is_duplicate and not read.is_secondary and not read.is_supplementary:
            if (read.pos > left_zero and read.pos < right_zero) or (read.next_reference_start > left_zero and read.next_reference_start < right_zero):
                if read.query_name not in matepairs:
                    matepairs[read.query_name] = read
                    
                else:
                    newname = str(uuid4())

                    mate = matepairs[read.query_name]

                    mate.query_name = newname
                    read.query_name = newname

                    if random.random() <= donor_norm_factor:
                        outbam.write(mate)
                        outbam.write(read)
                        paircount += 1

    outbam.close()

    logger.info('%s: using %d donor read pairs' % (mutid, paircount))

    return outbamfn


def fetch_read_names(args, chrom, start, end, svfrac=1.0):
    bamfile = pysam.AlignmentFile(args.bamFileName, 'rb')

    names = []

    for read in bamfile.fetch(chrom, start, end):
        if random.random() <= svfrac:
            names.append(read.query_name)

    return names


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

    if args.seed is not None: random.seed(int(args.seed) + int(bedline.split()[1]))

    mutid = '_'.join(map(str, bedline.split()[:4]))

    bamfile = pysam.Samfile(args.bamFileName, 'rb')
    reffile = pysam.Fastafile(args.refFasta)
    logfn = '_'.join(map(os.path.basename, bedline.split()[:4])) + ".log"
    logfile = open('addsv_logs_' + os.path.basename(args.outBamFile) + '/' + os.path.basename(args.outBamFile) + '_' + logfn, 'w')
    exclfile = args.tmpdir + '/' + '.'.join((mutid, 'exclude', str(uuid4()), 'txt'))
    exclude = open(exclfile, 'w')
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
                sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\t" + ' '.join(("copy number in sv region:",chrom,str(start),str(end),"=",str(cn))) + "\n")
                svfrac = svfrac/float(cn)
                assert svfrac <= 1.0, 'copy number from %s must be at least 1: %s' % (args.cnvfile, cnregion.strip())
                sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tadjusted default MAF: " + str(svfrac) + "\n")

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

    maxdfrac = 0.1 # FIXME make a parameter
    if dfrac > .1: 
        logger.warning("%s discordant fraction > %f aborting mutation!\n" % (mutid, maxdfrac))
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
                #print trn_refstart, trn_refend, user_trn_start, user_trn_end, int(args.maxlibsize)

            if trn_contig_end is None:
                logger.warning('%s contig does not cover user translocation end' % mutid)
                exact_success = False
                #print trn_refstart, trn_refend, user_trn_start, user_trn_end, int(args.maxlibsize)


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

            mutinfo[mutid] = "\t".join(('ins',chrom,str(refstart),str(refend),action,str(mutseq.length()),str(inspoint),str(insseqfile),str(tsdlen),str(svfrac)))
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

            mutinfo[mutid] = "\t".join(('trn',chrom,str(refstart),str(refend),action,str(trnpoint_1),trn_chrom,str(trn_refstart),str(trn_refend),str(trnpoint_2),str(svfrac)))
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
            raise ValueError("ERROR\t" + now() + "\t" + mutid + "\t: mutation not one of: INS,INV,DEL,DUP,TRN,BIGDEL,BIGINV,BIGDUP\n")

        logfile.write(">" + chrom + ":" + str(refstart) + "-" + str(refend) +" AFTER\n" + str(mutseq) + "\n")

    pemean, pesd = float(args.ismean), float(args.issd) 
    logger.info("%s set paired end mean distance: %f" % (mutid, pemean))
    logger.info("%s set paired end distance stddev: %f" % (mutid, pesd))


    # simulate reads
    (fq1, fq2) = runwgsim(maxcontig, mutseq.seq, svfrac, actions, exclude, pemean, pesd, args.tmpdir, err_rate=float(args.simerr), mutid=mutid, seed=args.seed, trn_contig=trn_maxcontig, rename=rename_reads)

    outreads = aligners.remap_fastq(args.aligner, fq1, fq2, args.refFasta, outbam_mutsfile, alignopts, mutid=mutid, threads=int(args.alignerthreads))

    if outreads == 0:
        logger.warning("%s outbam %s has no mapped reads!" % (mutid, outbam_mutsfile))
        return None, None, None

    logger.info("%s temporary bam: %s" % (mutid, outbam_mutsfile))

    exclude.close()
    bamfile.close()

    return outbam_mutsfile, exclfile, mutinfo




def main(args):
    logger.info("starting %s called with args: %s" % (sys.argv[0], ' '.join(sys.argv)))
    tmpbams = [] # temporary BAMs, each holds the realigned reads for one mutation
    exclfns = [] # 'exclude' files store reads to be removed from the original BAM due to deletions

    if not os.path.exists(args.bamFileName + '.bai'):
        logger.error("input bam must be indexed, not .bai file found for %s" % args.bamFileName)
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
        sys.stderr.write("\n")
        sys.exit(1)

    results = []
    pool = Pool(processes=int(args.procs))

    nmuts = 0

    if not os.path.exists(args.tmpdir):
        os.mkdir(args.tmpdir)
        logger.info("created tmp directory: %s" % args.tmpdir)

    if not os.path.exists('addsv_logs_' + os.path.basename(args.outBamFile)):
        os.mkdir('addsv_logs_' + os.path.basename(args.outBamFile))
        logger.info("created log directory: addsv_logs_%s" % os.path.basename(args.outBamFile))

    assert os.path.exists('addsv_logs_' + os.path.basename(args.outBamFile)), "could not create output directory!"
    assert os.path.exists(args.tmpdir), "could not create temporary directory!"

    bigdels = {}
    biginvs = {}
    bigdups = {}

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
                bdel_mutid = '_'.join(map(str, bedline.split()[:4]))
                bigdels[bdel_mutid] = (bdel_chrom, bdel_start, bdel_end, bdel_svfrac)

            # rewrite bigdup coords as translocation
            if mut_type == 'BIGDUP':
                bdup_svfrac = 1.0 # BIGDUP VAF is determined by donor bam

                if args.donorbam is None:
                    logger.warning('%s: using BIGDUP requires specifying a --donorbam and none was provided, using %s' % (bedline, args.bamFileName))
                    args.donorbam = args.bamFileName
                    continue

                bdup_chrom, bdup_start, bdup_end = bedline.split()[:3]
                bdup_start = int(bdup_start)
                bdup_end   = int(bdup_end)

                bdup_left_start = bdup_start
                bdup_left_end   = bdup_start

                bdup_right_start = bdup_end
                bdup_right_end   = bdup_end

                bedline = '%s %d %d BIGDUP %s %d %d %s %f' % (bdup_chrom, bdup_right_start, bdup_right_end, bdup_chrom, bdup_left_start, bdup_left_end, '++', bdup_svfrac)
                bdup_mutid = '_'.join(map(str, bedline.split()[:4]))
                bigdups[bdup_mutid] = (bdup_chrom, bdup_start, bdup_end, bdup_svfrac)

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
                result = pool.apply_async(makemut, [args, bedline, alignopts])
                results.append(result)

            else:
                for bedline in multi_part:
                    result = pool.apply_async(makemut, [args, bedline, alignopts])
                    results.append(result)

            nmuts += 1

    ## process the results of mutation jobs

    master_mutinfo = {}

    for result in results:
        tmpbam = None
        exclfn = None

        tmpbam, exclfn, mutinfo = result.get()

        if None not in (tmpbam, exclfn) and os.path.exists(tmpbam) and os.path.exists(exclfn):
            if bamreadcount(tmpbam) > 0:
                tmpbams.append(tmpbam)
                exclfns.append(exclfn)

                mutid = os.path.basename(tmpbam).split('.')[0]
                master_mutinfo[mutid] = mutinfo[mutid]

            else:
                os.remove(tmpbam)
                os.remove(exclfn)

    if len(tmpbams) == 0:
        logger.error("no succesful mutations")
        sys.exit()

    success_mutids = [os.path.basename(tmpbam).split('.')[0] for tmpbam in tmpbams]

    
    bigdel_excl = {}
    bigdup_add  = {}

    for mutid, mutinfo in master_mutinfo.items():
        # add additional excluded reads if bigdel(s) present
        if mutinfo.startswith('bigdel'):
            bdel_chrom, bdel_start, bdel_end, bdel_svfrac = bigdels[mutid]

            bdel_left_bnd = int(mutinfo.split()[3])
            bdel_right_bnd = int(mutinfo.split()[7])

            if bdel_left_bnd > bdel_right_bnd:
                bdel_left_bnd, bdel_right_bnd, bdel_right_bnd, bdel_left_bnd

            bigdel_excl[mutid] = fetch_read_names(args, bdel_chrom, bdel_left_bnd, bdel_right_bnd, svfrac=bdel_svfrac)

        if mutinfo.startswith('bigdup'):
            bdup_chrom, bdup_start, bdup_end, bdup_svfrac = bigdups[mutid]

            bdup_left_bnd = int((int(mutinfo.split()[7])+int(mutinfo.split()[8]))/2)
            bdup_right_bnd = int((int(mutinfo.split()[2])+int(mutinfo.split()[3]))/2)

            bigdup_add[mutid] = (bdup_chrom, bdup_left_bnd, bdup_right_bnd)


    biginv_pairs = dd(list)

    new_tmpbams = []

    for tmpbamfn in tmpbams:
        mutid = os.path.basename(tmpbamfn).split('.')[0]

        if mutid.endswith('BIGINV'):
            biginv_pairs[mutid].append(tmpbamfn)
        
        elif mutid.endswith('BIGDUP'):
            #print 'bigdup testing mutid:', mutid
            #print 'bigdup testing known mutids:', bigdup_add.keys()
            bdup_chrom, bdup_left_bnd, bdup_right_bnd = bigdup_add[mutid]

            bdup_left_bnd  = int(bdup_left_bnd)
            bdup_right_bnd = int(bdup_right_bnd)

            if bdup_left_bnd > bdup_right_bnd:
                bdup_left_bnd, bdup_right_bnd = bdup_right_bnd, bdup_left_bnd

            merged_bdup = add_donor_reads(args, mutid, tmpbamfn, bdup_chrom, bdup_left_bnd, bdup_right_bnd)

            new_tmpbams.append(merged_bdup) # TODO merge bigdup reads

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

            # add reads excluded due to BIGDEL if breakpoint was successful
            for bdel_mutid in bigdel_excl:
                if bdel_mutid in success_mutids:
                    for bdel_rn in bigdel_excl[bdel_mutid]:
                        exclout.write(bdel_rn+'\n')

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
        logger.info("cleaning up...")

        if not args.debug:
            if exclfn is not None:
                for exclfn in exclfns:
                    if os.path.isfile(exclfn):
                        os.remove(exclfn)

            for tmpbam in tmpbams:
                if os.path.isfile(tmpbam):
                    os.remove(tmpbam)
                if os.path.isfile(tmpbam + '.bai'):
                    os.remove(tmpbam + '.bai')

    else:
        if args.tagreads:
            from bamsurgeon.markreads import markreads
            tmp_tag_bam = 'tag.%s.bam' % str(uuid4())
            markreads(mergedtmp, tmp_tag_bam)
            move(tmp_tag_bam, mergedtmp)
            logger.info("tagged reads.")

        logger.info("swapping reads into original and writing to %s" % args.outBamFile)
        replace(args.bamFileName, mergedtmp, args.outBamFile, excl_merged, keepsecondary=args.keepsecondary, seed=args.seed)

        if not args.debug:
            os.remove(excl_merged)
            os.remove(mergedtmp)

            for exclfn in exclfns:
                if os.path.isfile(exclfn):
                    os.remove(exclfn)

            for tmpbam in tmpbams:
                if os.path.isfile(tmpbam):
                    os.remove(tmpbam)
                if os.path.isfile(tmpbam + '.bai'):
                    os.remove(tmpbam + '.bai')

        logger.info("done.")

    var_basename = '.'.join(os.path.basename(args.varFileName).split('.')[:-1])
    bam_basename = '.'.join(os.path.basename(args.outBamFile).split('.')[:-1])

    vcf_fn = bam_basename + '.addindel.' + var_basename + '.vcf'

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
    parser.add_argument('--seed', default=None,
                        help='seed random number generation')
    parser.add_argument('--allowN', action='store_true', default=False,
                        help='allow N in contigs, replace with A and warn user (default: drop mutation)')
    args = parser.parse_args()
    main(args)

