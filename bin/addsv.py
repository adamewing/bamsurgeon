#!/usr/bin/env python

#from __future__ import print_function

import re
import os
import sys
import random
import subprocess
import traceback
import argparse
import pysam
import bamsurgeon.replacereads as rr
import bamsurgeon.asmregion as ar
import bamsurgeon.mutableseq as ms
import bamsurgeon.aligners as aligners

from bamsurgeon.common import *
from uuid import uuid4
from time import sleep
from shutil import move
from math import sqrt
from itertools import izip
from collections import Counter
from multiprocessing import Pool

import logging
FORMAT = '%(levelname)s %(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
sys.stderr = os.fdopen(sys.stderr.fileno(), 'w', 0)


def runwgsim(contig, newseq, svfrac, svtype, exclude, pemean, pesd, tmpdir, mutid='null', err_rate=0.0, seed=None, trn_contig=None):
    ''' wrapper function for wgsim
    '''

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

    refstart = start + tgtstart
    refend = start + tgtend

    if refstart > refend:
        refstart, refend = refend, refstart

    return contig, refseq, alignstats, refstart, refend, qrystart, qryend, tgtstart, tgtend


def fetch_read_names(args, chrom, start, end, svfrac=1.0):
    bamfile = pysam.Samfile(args.bamFileName, 'rb')

    names = []

    for read in bamfile.fetch(chrom, start, end):
        if random.random() <= svfrac:
            names.append(read.qname)

    return names


def makemut(args, bedline, alignopts):

    if args.seed is not None: random.seed(int(args.seed) + int(bedline.strip().split()[1]))

    mutid = '_'.join(map(str, bedline.strip().split()[:4]))

    try:
        bamfile = pysam.Samfile(args.bamFileName, 'rb')
        reffile = pysam.Fastafile(args.refFasta)
        logfn = '_'.join(map(os.path.basename, bedline.strip().split()[:4])) + ".log"
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

        c = bedline.strip().split()
        chrom  = c[0]
        start  = int(c[1])
        end    = int(c[2])
        araw   = c[3:] # INV, DEL, INS,  DUP, TRN
 
        # translocation specific
        trn_chrom = None
        trn_start = None
        trn_end   = None

        is_transloc = c[3] in ('TRN', 'BIGDEL')

        if is_transloc:
            araw = [c[3]]
            if len(c) > 7:
                araw += c[7:]

            start -= 3000
            end   += 3000
            if start < 0: start = 0

            trn_chrom = c[4]
            trn_start = int(c[5]) - 3000
            trn_end   = int(c[5]) + 3000
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
                    assert svfrac <= 1.0, 'copy number from %s must be at least 1: %s' % (args.cnvfile, snregion.stri[()])
                    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tadjusted default MAF: " + str(svfrac) + "\n")

        logger.info("%s interval: %s" % (mutid, bedline.strip()))
        logger.info("%s length: %d" % (mutid, (end-start)))

       # modify start and end if interval is too short
        minctglen = int(args.minctglen)

        # adjust if minctglen is too short
        if minctglen < 3*int(args.maxlibsize):
            minctglen = 3*int(args.maxlibsize)

        if end-start < minctglen:
            adj   = minctglen - (end-start)
            start = start - adj/2
            end   = end + adj/2

            logger.info("%s note: interval size was too short, adjusted: %s:%d-%d" % (mutid, chrom,start,end))

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

        if is_transloc: trn_mutseq = ms.MutableSeq(trn_maxcontig.seq)

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
                if len(a) > 1:
                    dsize = float(a[1])
                    if dsize > 1.0: # if DEL size is not a fraction, interpret as bp
                        # since DEL 1 is default, if DEL 1 is specified, interpret as 1 bp deletion
                        dlen = int(dsize)
                        dsize = 1.0
                else:
                    dsize = 1.0

                if len(a) > 2: # VAF
                    svfrac = float(a[2])/cn

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

            if action == 'INS':
                inspoint = mutseq.length()/2
                if ins_motif is not None:
                    inspoint = mutseq.find_site(ins_motif, left_trim=int(args.maxlibsize), right_trim=int(args.maxlibsize))

                    if inspoint < int(args.maxlibsize) or inspoint > mutseq.length() - int(args.maxlibsize):
                        logger.info("%s picked midpoint, no cutsite found" % mutid)
                        inspoint = mutseq.length()/2

                if insseqfile: # seq in file
                    if insseqfile == 'RND':
                        assert args.inslib is not None # insertion library needs to exist
                        insseqfile = random.choice(args.inslib.keys())
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
                mutseq.inversion(invstart,invend)

                mutinfo[mutid] = "\t".join(('inv',chrom,str(refstart),str(refend),action,str(mutseq.length()),str(invstart),str(invend),str(svfrac)))
                logfile.write(mutinfo[mutid] + "\n")

            elif action == 'DEL':
                delstart = int(args.maxlibsize)
                delend = mutseq.length() - delstart
                if dlen == 0: # bp size not specified, delete fraction of contig
                    dlen = int((float(delend-delstart) * dsize)+0.5) 

                dadj = delend-delstart-dlen
                if dadj < 0:
                    dadj = 0
                    logger.warning("%s warning: deletion of length 0" % mutid)

                delstart += dadj/2
                delend   -= dadj/2

                mutseq.deletion(delstart,delend)

                mutinfo[mutid] = "\t".join(('del',chrom,str(refstart),str(refend),action,str(mutseq.length()),str(delstart),str(delend),str(dlen),str(svfrac)))
                logfile.write(mutinfo[mutid] + "\n")

            elif action == 'DUP':
                dupstart = int(args.maxlibsize)
                dupend = mutseq.length() - dupstart
                mutseq.duplication(dupstart,dupend,ndups)

                mutinfo[mutid] = "\t".join(('dup',chrom,str(refstart),str(refend),action,str(mutseq.length()),str(dupstart),str(dupend),str(ndups),str(svfrac)))
                logfile.write(mutinfo[mutid] + "\n")

            elif action == 'TRN':
                mutseq.fusion(mutseq.length()/2, trn_mutseq, trn_mutseq.length()/2, flip1=trn_left_flip, flip2=trn_right_flip)

                mutinfo[mutid] = "\t".join(('trn',chrom,str(refstart),str(refend),action,str(mutseq.length()),trn_chrom,str(trn_refstart),str(trn_refend),str(trn_mutseq.length()),str(svfrac)))
                logfile.write(mutinfo[mutid] + "\n")

            elif action == 'BIGDEL':
                mutseq.fusion(mutseq.length()/2, trn_mutseq, trn_mutseq.length()/2)

                mutinfo[mutid] = "\t".join(('bigdel',chrom,str(refstart),str(refend),action,str(mutseq.length()),trn_chrom,str(trn_refstart),str(trn_refend),str(trn_mutseq.length()),str(svfrac)))
                logfile.write(mutinfo[mutid] + "\n")

            elif action == 'BIGDUP':
                raise NotImplementedError()

            elif action == 'BIGINV':
                raise NotImplementedError()

            else:
                raise ValueError("ERROR\t" + now() + "\t" + mutid + "\t: mutation not one of: INS,INV,DEL,DUP,TRN,BIGDEL\n")

            logfile.write(">" + chrom + ":" + str(refstart) + "-" + str(refend) +" AFTER\n" + str(mutseq) + "\n")

        pemean, pesd = float(args.ismean), float(args.issd) 
        logger.info("%s set paired end mean distance: %f" % (mutid, pemean))
        logger.info("%s set paired end distance stddev: %f" % (mutid, pesd))

        # simulate reads
        (fq1, fq2) = runwgsim(maxcontig, mutseq.seq, svfrac, actions, exclude, pemean, pesd, args.tmpdir, err_rate=float(args.simerr), mutid=mutid, seed=args.seed, trn_contig=trn_maxcontig)

        outreads = aligners.remap_fastq(args.aligner, fq1, fq2, args.refFasta, outbam_mutsfile, alignopts, mutid=mutid, threads=1)

        if outreads == 0:
            logger.warning("%s outbam %s has no mapped reads!" % (mutid, outbam_mutsfile))
            return None, None, None

        logger.info("%s temporary bam: %s" % (mutid, outbam_mutsfile))

        exclude.close()
        bamfile.close()

        return outbam_mutsfile, exclfile, mutinfo

    except Exception, e:
        sys.stderr.write("*"*60 + "\nencountered error in mutation spikein: " + bedline + "\n")
        traceback.print_exc(file=sys.stderr)
        sys.stderr.write("*"*60 + "\n")
        return None, None, None


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
    except Exception, e:
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

    with open(args.varFileName, 'r') as varfile:
        for bedline in varfile:
            if re.search('^#',bedline):
                continue

            if args.maxmuts and nmuts >= int(args.maxmuts):
                break

            # rewrite bigdel coords as translocation
            if bedline.strip().split()[3] == 'BIGDEL':
                bd_svfrac=1.0
                if len(bedline.strip().split()) == 5:
                    bd_svfrac = float(bedline.strip().split()[-1])

                # rewrite bigdel coords as translocation
                bd_chrom, bd_start, bd_end = bedline.strip().split()[:3]
                bd_start = int(bd_start)
                bd_end = int(bd_end)

                bd_left_start = bd_start - 2000
                bd_left_end = bd_start + 2000

                bd_right_start = bd_end - 2000
                bd_right_end = bd_end + 2000

                bedline = '%s %d %d BIGDEL %s %d %d %s %f' % (bd_chrom, bd_left_start, bd_left_end, bd_chrom, bd_right_start, bd_right_end, '++', bd_svfrac)
                bd_mutid = '_'.join(map(str, bedline.strip().split()[:4]))
                bigdels[bd_mutid] = (bd_chrom, bd_start, bd_end, bd_svfrac)

            # submit each mutation as its own thread
            result = pool.apply_async(makemut, [args, bedline, alignopts])
            results.append(result)

            nmuts += 1
            if args.delay is not None:
                sleep(int(args.delay))

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

    # add additional excluded reads if bigdel(s) present
    bigdel_excl = {}
    for mutid, mutinfo in master_mutinfo.iteritems():
        if mutinfo.startswith('bigdel'):
            bd_chrom, bd_start, bd_end, bd_svfrac = bigdels[mutid]

            bd_left_bnd = int(mutinfo.split()[3])
            bd_right_bnd = int(mutinfo.split()[7])

            bigdel_excl[mutid] = fetch_read_names(args, bd_chrom, bd_left_bnd, bd_right_bnd, svfrac=bd_svfrac)

    logger.info("tmpbams: %s" % tmpbams)
    logger.info("exclude: %s" % exclfns)

    excl_merged = 'addsv.exclude.final.' + str(uuid4()) + '.txt'
    mergedtmp = 'addsv.mergetmp.final.' + str(uuid4()) + '.bam'

    logger.info("merging exclude files into %s" % excl_merged)
    exclout = open(excl_merged, 'w')
    for exclfn in exclfns:
        with open(exclfn, 'r') as excl:
            for line in excl:
                exclout.write(line)

            # add reads excluded due to BIGDEL if breakpoint was successful
            for bd_mutid in bigdel_excl:
                if bd_mutid in success_mutids:
                    for bd_rn in bigdel_excl[bd_mutid]:
                        exclout.write(bd_rn+'\n')

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
    parser.add_argument('--minctglen', dest='minctglen', default=3000,
                        help="pad input intervals out to a minimum length for contig generation (default=3000)")
    parser.add_argument('-n', dest='maxmuts', default=None,
                        help="maximum number of mutations to make")
    parser.add_argument('-c', '--cnvfile', dest='cnvfile', default=None, 
                        help="tabix-indexed list of genome-wide absolute copy number values (e.g. 2 alleles = no change)")
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
    parser.add_argument('--delay', default=None, 
                        help='time delay between jobs (try to avoid thrashing disks)')
    parser.add_argument('--noref', action='store_true', default=False, 
                        help="do not perform reference based assembly")
    parser.add_argument('--aligner', default='backtrack',
                        help='supported aligners: ' + ','.join(aligners.supported_aligners_fastq))
    parser.add_argument('--alignopts', default=None,
                        help='aligner-specific options as comma delimited list of option1:value1,option2:value2,...')
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

