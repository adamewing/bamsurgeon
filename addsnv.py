#!/usr/bin/env python

import sys
import pysam
import argparse
import random
import subprocess
import os
import bs.replacereads as rr
import traceback
import bs.aligners as aligners

from bs.common import *
from uuid import uuid4
from shutil import move
from re import sub
from multiprocessing import Pool, Value, Lock


sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
sys.stderr = os.fdopen(sys.stderr.fileno(), 'w', 0)


def mut(base, altbase):
    """ change base to something different
    """

    bases = ('A','T','C','G')
    base = base.upper()
    if base not in bases or (altbase is not None and altbase not in ['A','T','G','C']):
        raise ValueError("ERROR\t" + now() + "\tbase passed to mut(): " + str(base) + " not one of (A,T,C,G)\n")

    if altbase is not None:
        return altbase

    else:
        alt = base
        while alt == base:
            alt = bases[int(random.uniform(0,4))]
        return alt

def countReadCoverage(bam,chrom,start,end):
    """ calculate coverage of aligned reads over region
    """

    coverage = []
    start = int(start)
    end = int(end)
    for i in range(end-start+1):
        coverage.append(0.0)

    i = 0
    if chrom in bam.references:
        for pcol in bam.pileup(chrom,start,end):
            n = 0
            if pcol.pos >= start and pcol.pos <= end:
                for read in pcol.pileups:
                    if read.alignment.mapq >= 0 and not read.alignment.is_duplicate:
                        n += 1
                coverage[i] = n
                i += 1

    return coverage


def countBaseAtPos(bamfile, chrom, pos, mutid='null'):
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
            sys.stderr.write("WARN\t" + now() + "\t" + mutid + "\tmpileup failed, no coverage for base: " + chrom + ":" + str(pos))
            return []
    bases = []
    if pileup:
        for b in pileup:
            if b in ['A','T','C','G']:
                bases.append(b)

    return bases


def replace(origbamfile, mutbamfile, outbamfile):
    ''' open .bam file and call replacereads
    '''
    origbam = pysam.Samfile(origbamfile, 'rb')
    mutbam  = pysam.Samfile(mutbamfile, 'rb')
    outbam  = pysam.Samfile(outbamfile, 'wb', template=origbam)

    rr.replaceReads(origbam, mutbam, outbam, keepqual=True)

    origbam.close()
    mutbam.close()
    outbam.close()



def makemut(args, chrom, start, end, vaf, altbase, avoid, alignopts):
    mutid = chrom + '_' + str(start) + '_' + str(end) + '_' + str(vaf) + '_' + str(altbase)
    try:
        bamfile = pysam.Samfile(args.bamFileName, 'rb')
        bammate = pysam.Samfile(args.bamFileName, 'rb') # use for mates to avoid iterator problems
        reffile = pysam.Fastafile(args.refFasta)
        tmpbams = []

        snvfrac = float(args.snvfrac)

        mutpos = int(random.uniform(start,end+1)) # position of mutation in genome
        refbase = reffile.fetch(chrom,mutpos-1,mutpos)

        if altbase == refbase.upper():
            sys.stderr.write("WARN\t" + now() + "\t" + mutid + "\tspecified ALT base matches reference, skipping mutation\n")
            return None

        try:
            mutbase = mut(refbase, altbase)
        except ValueError as e:
            sys.stderr.write("WARN\t" + now() + "\t" + mutid + "\t" + ' '.join(("skipped site:",chrom,str(start),str(end),"due to N base:",str(e),"\n")))
            return None

        mutstr = refbase + "-->" + mutbase

        # optional CNV file
        cnv = None
        if (args.cnvfile):
            cnv = pysam.Tabixfile(args.cnvfile, 'r')

        log = open('addsnv_logs_' + os.path.basename(args.outBamFile) + '/' + os.path.basename(args.outBamFile) + "." + "_".join((chrom,str(start),str(end))) + ".log",'w')

        # keep a list of reads to modify - use hash to keep unique since each
        # read will be visited as many times as it has bases covering the region
        outreads = {}
        mutreads = {} # same keys as outreads
        mutmates = {} # same keys as outreads, keep track of mates
        numunmap = 0
        hasSNP = False
        tmpoutbamname = args.tmpdir + "/" + mutid + ".tmpbam." + str(uuid4()) + ".bam"
        print "INFO\t" + now() + "\t" + mutid + "\tcreating tmp bam: ",tmpoutbamname
        outbam_muts = pysam.Samfile(tmpoutbamname, 'wb', template=bamfile)
        maxfrac = 0.0

        for pcol in bamfile.pileup(reference=chrom,start=mutpos,end=mutpos+1):
            # this will include all positions covered by a read that covers the region of interest
            if pcol.pos: #> start and pcol.pos <= end:
                refbase = reffile.fetch(chrom,pcol.pos-1,pcol.pos)
                basepile = ''
                for pread in pcol.pileups:
                    if avoid is not None and pread.alignment.qname in avoid:
                        print "WARN\t" + now() + "\t" + mutid + "\tdropped mutation due to read in --avoidlist", pread.alignment.qname
                        os.remove(tmpoutbamname)
                        return None

                    if not pread.alignment.is_secondary and bin(pread.alignment.flag & 2048) != bin(2048): # only consider primary alignments
                        basepile += pread.alignment.seq[pread.query_position-1]
                        pairname = 'F' # read is first in pair
                        if pread.alignment.is_read2:
                            pairname = 'S' # read is second in pair
                        if not pread.alignment.is_paired:
                            pairname = 'U' # read is unpaired

                        extqname = ','.join((pread.alignment.qname,str(pread.alignment.pos),pairname))

                        if pcol.pos == mutpos:
                            if not pread.alignment.is_secondary and bin(pread.alignment.flag & 2048) != bin(2048) and not pread.alignment.mate_is_unmapped:
                                outreads[extqname] = pread.alignment
                                mutbases = list(pread.alignment.seq)
                                mutbases[pread.query_position-1] = mutbase
                                mutread = ''.join(mutbases)
                                mutreads[extqname] = mutread
                                mate = None
                                if not args.single:
                                    try:
                                        mate = bammate.mate(pread.alignment)
                                    except:
                                        print "WARN\t" + now() + "\t" + mutid + "\twarning: no mate for", pread.alignment.qname
                                        if args.requirepaired:
                                            print "WARN\t" + now() + "\t" + mutid + "\tskipped mutation due to --requirepaired"
                                            os.remove(tmpoutbamname)
                                            return None
                                mutmates[extqname] = mate
                                log.write(" ".join(('read',extqname,mutread,"\n")))
                            else:
                                numunmap += 1

                            if len(mutreads) > int(args.maxdepth):
                                sys.stderr.write("WARN\t" + now() + "\t" + mutid + "\tdepth at site is greater than cutoff, aborting mutation.\n")
                                outbam_muts.close()
                                os.remove(tmpoutbamname)
                                return None

                # make sure region doesn't have any changes that are likely SNPs
                # (trying to avoid messing with haplotypes)
                    
                basepile = countBaseAtPos(args.bamFileName,chrom,pcol.pos,mutid=mutid)
                if basepile:
                    majb = majorbase(basepile)
                    minb = minorbase(basepile)

                    frac = float(minb[1])/(float(majb[1])+float(minb[1]))
                    if minb[0] == majb[0]:
                        frac = 0.0
                    if frac > maxfrac:
                        maxfrac = frac
                    if (not args.ignoresnps) and (frac > snvfrac or majb[0].upper() != refbase.upper()):
                        print "WARN\t" + now() + "\t" + mutid + "\tdropped for proximity to SNP, nearby SNP MAF:",frac,"maxfrac:",snvfrac
                        hasSNP = True
                else:
                    print "WARN\t" + now() + "\t" + mutid + "\tcould not pileup for region:",chrom,pcol.pos
                    if not args.ignorepileup:
                        hasSNP = True

        # pick reads to change
        readlist = []
        for extqname,read in outreads.iteritems():
            if read.seq != mutreads[extqname]:
                readlist.append(extqname)

        print "INFO\t" + now() + "\t" + mutid + "\tlen(readlist): " + str(len(readlist))
        random.shuffle(readlist)

        if len(readlist) < int(args.mindepth):
            print "WARN\t" + now() + "\t" + mutid + "\ttoo few reads in region (" + str(len(readlist)) + ") skipping..."
            outbam_muts.close()
            os.remove(tmpoutbamname)
            return None

        if vaf is None:
            vaf = float(args.mutfrac) # default minor allele freq if not otherwise specified
        if cnv: # cnv file is present
            if chrom in cnv.contigs:
                for cnregion in cnv.fetch(chrom,start,end):
                    cn = float(cnregion.strip().split()[3]) # expect chrom,start,end,CN
                    print "INFO\t" + now() + "\t" + mutid + "\t" + ' '.join(("copy number in snp region:",chrom,str(start),str(end),"=",str(cn))) + "\n"
                    if float(cn) > 0.0:
                        vaf = 1.0/float(cn)
                    else:
                        vaf = 0.0
                    print "adjusted VAF: " + str(vaf) + "\n"
        else:
            print "INFO\t" + now() + "\t" + mutid + "\tselected VAF: " + str(vaf) + "\n"

        ######
        lastread = int(len(readlist)*vaf)

        # pick at least args.minmutreads if possible
        if lastread < int(args.minmutreads):
            if len(readlist) > int(args.minmutreads):
                lastread = int(args.minmutreads)
                sys.stdout.write("WARN\t" + now() + "\t" + mutid + "\tforced " + str(lastread) + " reads.\n")
            else:
                print "WARN\t" + now() + "\t" + mutid + "\tdropped site with fewer reads than --minmutreads"
                os.remove(tmpoutbamname)
                return None

        readlist = readlist[0:lastread] 

        ######

        print "INFO\t" + now() + "\t" + mutid + "\tpicked:",str(len(readlist))

        wrote = 0
        nmut = 0
        # change reads from .bam to mutated sequences
        for extqname,read in outreads.iteritems():
            if read.seq != mutreads[extqname]:
                if not args.nomut and extqname in readlist:
                    qual = read.qual # changing seq resets qual (see pysam API docs)
                    read.seq = mutreads[extqname] # make mutation
                    read.qual = qual
                    nmut += 1
            if not hasSNP or args.force:
                wrote += 1
                outbam_muts.write(read)
                if args.single:
                    outbam_muts.write(mutmates[extqname])
                elif mutmates[extqname] is not None:
                    outbam_muts.write(mutmates[extqname])
        print "INFO\t" + now() + "\t" + mutid + "\twrote: ",wrote,"mutated:",nmut

        if not hasSNP or args.force:
            outbam_muts.close()

            aligners.remap_bam(args.aligner, tmpoutbamname, args.refFasta, alignopts, mutid=mutid, paired=(not args.single), samtofastq=args.samtofastq)

            outbam_muts = pysam.Samfile(tmpoutbamname,'rb')
            coverwindow = 1
            incover  = countReadCoverage(bamfile,chrom,mutpos-coverwindow,mutpos+coverwindow)
            outcover = countReadCoverage(outbam_muts,chrom,mutpos-coverwindow,mutpos+coverwindow)

            avgincover  = float(sum(incover))/float(len(incover)) 
            avgoutcover = float(sum(outcover))/float(len(outcover))

            print "INFO\t" + now() + "\t" + mutid + "\tavgincover: " + str(avgincover) + " avgoutcover: " + str(avgoutcover)

            spikein_snvfrac = 0.0
            if wrote > 0:
                spikein_snvfrac = float(nmut)/float(wrote)

            # qc cutoff for final snv depth 
            if (avgoutcover > 0 and avgincover > 0 and avgoutcover/avgincover >= float(args.coverdiff)) or args.force:
                tmpbams.append(tmpoutbamname)
                snvstr = chrom + ":" + str(start) + "-" + str(end) + " (VAF=" + str(vaf) + ")"
                log.write("\t".join(("snv",snvstr,str(mutpos),mutstr,str(avgoutcover),str(avgoutcover),str(spikein_snvfrac),str(maxfrac)))+"\n")
            else:
                outbam_muts.close()
                os.remove(tmpoutbamname)
                if os.path.exists(tmpoutbamname + '.bai'):
                    os.remove(tmpoutbamname + '.bai')

                return None

        outbam_muts.close()
        bamfile.close()
        bammate.close()
        log.close() 

        return tmpbams

    except Exception, e:
        sys.stderr.write("*"*60 + "\nERROR\t" + now() + "\tencountered error in mutation spikein: " + mutid + "\n")
        traceback.print_exc(file=sys.stdout)
        sys.stderr.write("*"*60 + "\n")
        if os.path.exists(tmpoutbamname):
            os.remove(tmpoutbamname)
        if os.path.exists(tmpoutbamname + '.bai'):
            os.remove(tmpoutbamname + '.bai')
        return None

def main(args):
    print "INFO\t" + now() + "\tstarting " + sys.argv[0] + " called with args: " + ' '.join(sys.argv) + "\n"
    bedfile = open(args.varFileName, 'r')
    reffile = pysam.Fastafile(args.refFasta)

    if not os.path.exists(args.bamFileName + '.bai'):
        sys.stderr.write("ERROR\t" + now() + "\tinput bam must be indexed, not .bai file found for " + args.bamFileName + " \n")
        sys.exit(1)

    alignopts = {}
    if args.alignopts is not None:
        alignopts = dict([o.split(':') for o in args.alignopts.split(',')])

    aligners.checkoptions(args.aligner, alignopts, args.samtofastq)

    # load readlist to avoid, if specified
    avoid = None
    if args.avoidreads is not None:
        avoid = dictlist(args.avoidreads)

    # make a temporary file to hold mutated reads
    outbam_mutsfile = "addsnv." + str(uuid4()) + ".muts.bam"
    bamfile = pysam.Samfile(args.bamFileName, 'rb')
    outbam_muts = pysam.Samfile(outbam_mutsfile, 'wb', template=bamfile)
    outbam_muts.close()
    bamfile.close()
    tmpbams = []

    if not os.path.exists(args.tmpdir):
        os.mkdir(args.tmpdir)
        print "INFO\t" + now() + "\tcreated tmp directory: " + args.tmpdir

    if not os.path.exists('addsnv_logs_' + os.path.basename(args.outBamFile)):
        os.mkdir('addsnv_logs_' + os.path.basename(args.outBamFile))
        print "INFO\t" + now() + "\tcreated directory: addsnv_logs_" + os.path.basename(args.outBamFile)

    assert os.path.exists('addsnv_logs_' + os.path.basename(args.outBamFile)), "could not create output directory!"
    assert os.path.exists(args.tmpdir), "could not create temporary directory!"

    pool = Pool(processes=int(args.procs))
    results = []

    ntried = 0
    for bedline in bedfile:
        if ntried < int(args.numsnvs) or int(args.numsnvs) == 0:
            c = bedline.strip().split()
            chrom   = c[0]
            start   = int(c[1])
            end     = int(c[2])
            vaf     = None
            altbase = None

            # VAF is 4th column, if present
            if len(c) > 3:
                vaf = float(c[3])

            # ALT is 5th column, if present
            if len(c) == 5:
                altbase = c[4].upper()
                assert altbase in ['A','T','C','G'], "ERROR:\t" + now() + "\tALT " + altbase + " not A, T, C, or G!\n"

            # make mutation (submit job to thread pool)
            result = pool.apply_async(makemut, [args, chrom, start, end, vaf, altbase, avoid, alignopts])
            results.append(result)
            ntried += 1

    for result in results:
        tmpbamlist = result.get()
        if tmpbamlist is not None:
            for tmpbam in tmpbamlist:
                tmpbams.append(tmpbam)

    # merge tmp bams
    if len(tmpbams) == 1:
        move(tmpbams[0],outbam_mutsfile)
    elif len(tmpbams) > 1:
        mergebams(tmpbams,outbam_mutsfile,maxopen=int(args.maxopen))

    bedfile.close()

    # cleanup
    for bam in tmpbams:
        if os.path.exists(bam):
            os.remove(bam)
        if os.path.exists(bam + '.bai'):
            os.remove(bam + '.bai')

    if args.skipmerge:
        print "INFO\t" + now() + "\tskipping merge, plase merge reads from", outbam_mutsfile, "manually."
    else:
        print "INFO\t" + now() + "\tdone making mutations, merging mutations into", args.bamFileName, "-->", args.outBamFile
        replace(args.bamFileName, outbam_mutsfile, args.outBamFile)

        #cleanup
        os.remove(outbam_mutsfile)


def run():
    parser = argparse.ArgumentParser(description='adds SNVs to reads, outputs modified reads as .bam along with mates')
    parser.add_argument('-v', '--varfile', dest='varFileName', required=True, help='Target regions to try and add a SNV, as BED')
    parser.add_argument('-f', '--bamfile', dest='bamFileName', required=True, help='sam/bam file from which to obtain reads')
    parser.add_argument('-r', '--reference', dest='refFasta', required=True, help='reference genome, fasta indexed with bwa index -a stdsw _and_ samtools faidx')
    parser.add_argument('-o', '--outbam', dest='outBamFile', required=True, help='.bam file name for output')
    parser.add_argument('-s', '--snvfrac', dest='snvfrac', default=1, help='maximum allowable linked SNP MAF (for avoiding haplotypes) (default = 1)')
    parser.add_argument('-m', '--mutfrac', dest='mutfrac', default=0.5, help='allelic fraction at which to make SNVs (default = 0.5)')
    parser.add_argument('-n', '--numsnvs', dest='numsnvs', default=0, help="maximum number of mutations to try (default: entire input)")
    parser.add_argument('-c', '--cnvfile', dest='cnvfile', default=None, help="tabix-indexed list of genome-wide absolute copy number values (e.g. 2 alleles = no change)")
    parser.add_argument('-d', '--coverdiff', dest='coverdiff', default=0.9, help="allow difference in input and output coverage (default=0.9)")
    parser.add_argument('-p', '--procs', dest='procs', default=1, help="split into multiple processes (default=1)")
    parser.add_argument('--samtofastq', default=None, help='path to picard SamToFastq.jar, required for most aligners')
    parser.add_argument('--mindepth', default=10, help='minimum read depth to make mutation (default = 10)')
    parser.add_argument('--maxdepth', default=2000, help='maximum read depth to make mutation (default = 2000)')
    parser.add_argument('--minmutreads', default=3, help='minimum number of mutated reads to output per site')
    parser.add_argument('--avoidreads', default=None, help='file of read names to avoid (mutations will be skipped if overlap)')
    parser.add_argument('--nomut', action='store_true', default=False, help="dry run")
    parser.add_argument('--ignoresnps', action='store_true', default=False, help="make mutations even if there are non-reference alleles sharing the relevant reads")
    parser.add_argument('--force', action='store_true', default=False, help="force mutation to happen regardless of nearby SNP or low coverage")
    parser.add_argument('--single', action='store_true', default=False, help="input BAM is simgle-ended (default is paired-end)")
    parser.add_argument('--maxopen', dest='maxopen', default=1000, help="maximum number of open files during merge (default 1000)")
    parser.add_argument('--requirepaired', action='store_true', default=False, help='skip mutations if unpaired reads are present')
    parser.add_argument('--skipmerge', action='store_true', default=False, help="final output is tmp file to be merged")
    parser.add_argument('--ignorepileup', action='store_true', default=False, help="do not check pileup depth in mutation regions")
    parser.add_argument('--aligner', default='backtrack', help='supported aligners: ' + ','.join(aligners.supported_aligners_bam))
    parser.add_argument('--alignopts', default=None, help='aligner-specific options as comma delimited list of option1:value1,option2:value2,...')
    parser.add_argument('--tmpdir', default='addsnv.tmp', help='temporary directory (default=addsnv.tmp)')
    args = parser.parse_args()
    main(args)

if __name__ == '__main__':
    run()
