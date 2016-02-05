#!/usr/bin/env python

import sys
import pysam
import argparse
import random
import subprocess
import os
import bamsurgeon.replacereads as rr
import bamsurgeon.aligners as aligners
import bamsurgeon.mutation as mutation
import traceback

from bamsurgeon.common import *
from uuid import uuid4
from re import sub
from shutil import move
from multiprocessing import Pool
from collections import defaultdict as dd

sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
sys.stderr = os.fdopen(sys.stderr.fileno(), 'w', 0)


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


def replace(origbamfile, mutbamfile, outbamfile, seed=None):
    ''' open .bam file and call replacereads
    '''
    origbam = pysam.Samfile(origbamfile, 'rb')
    mutbam  = pysam.Samfile(mutbamfile, 'rb')
    outbam  = pysam.Samfile(outbamfile, 'wb', template=origbam)

    rr.replaceReads(origbam, mutbam, outbam, keepqual=True, seed=seed)

    origbam.close()
    mutbam.close()
    outbam.close()


def get_mutstr(chrom, start, end, ins, ref):
    return 'FIX get_mutstr'


def dictlist(fn):
    d = {}
    with open(fn, 'r') as inlist:
        for name in inlist:
            d[name.strip()] = True
    return d


def makemut(args, chrom, start, end, vaf, ins, avoid, alignopts):
    ''' is ins is a sequence, it will is inserted at start, otherwise delete from start to end'''

    if args.seed is not None: random.seed(int(args.seed) + int(start))

    mutid = chrom + '_' + str(start) + '_' + str(end) + '_' + str(vaf)
    if ins is None:
        mutid += ':DEL'
    else:
        mutid += ':INS:' + ins

    try:
        bamfile = pysam.Samfile(args.bamFileName, 'rb')
        bammate = pysam.Samfile(args.bamFileName, 'rb') # use for mates to avoid iterator problems
        reffile = pysam.Fastafile(args.refFasta)
        tmpbams = []

        is_insertion = ins is not None
        is_deletion  = ins is None

        snvfrac = float(args.snvfrac)

        mutstr = get_mutstr(chrom, start, end, ins, reffile)

        del_ln = 0
        if is_deletion:
            del_ln = end-start

        mutpos = start
        mutpos_list = [start]

        # optional CNV file
        cnv = None
        if (args.cnvfile):
            cnv = pysam.Tabixfile(args.cnvfile, 'r')

        log = open('addindel_logs_' + os.path.basename(args.outBamFile) + '/' + os.path.basename(args.outBamFile) + "." + "_".join((chrom,str(start),str(end))) + ".log",'w')

        tmpoutbamname = args.tmpdir + "/" + mutid + ".tmpbam." + str(uuid4()) + ".bam"
        print "INFO\t" + now() + "\t" + mutid + "\tcreating tmp bam: ",tmpoutbamname #DEBUG
        outbam_muts = pysam.Samfile(tmpoutbamname, 'wb', template=bamfile)

        mutfail, hasSNP, maxfrac, outreads, mutreads, mutmates = mutation.mutate(args, log, bamfile, bammate, chrom, mutpos, mutpos+del_ln+1, mutpos_list, avoid=avoid, mutid_list=[mutid], is_insertion=is_insertion, is_deletion=is_deletion, ins_seq=ins, reffile=reffile, indel_start=start, indel_end=end)

        if mutfail:
            outbam_muts.close()
            os.remove(tmpoutbamname)
            return None

        # pick reads to change
        readlist = []
        for extqname,read in outreads.iteritems():
            if read.seq != mutreads[extqname]:
                readlist.append(extqname)

        print "len(readlist):",str(len(readlist))
        readlist.sort()
        random.shuffle(readlist)

        if len(readlist) < int(args.mindepth):
            sys.stderr.write("WARN\t" + now() + "\t" + mutid + "\tskipped, too few reads in region: " + str(len(readlist)) + "\n")
            outbam_muts.close()
            os.remove(tmpoutbamname)
            return None

        if vaf is None:
            vaf = float(args.mutfrac) # default minor allele freq if not otherwise specified
        if cnv: # cnv file is present
            if chrom in cnv.contigs:
                for cnregion in cnv.fetch(chrom,start,end):
                    cn = float(cnregion.strip().split()[3]) # expect chrom,start,end,CN
                    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\t" + ' '.join(("copy number in snp region:",chrom,str(start),str(end),"=",str(cn))) + "\n")
                    if float(cn) > 0.0:
                        vaf = 1.0/float(cn)
                    else:
                        vaf = 0.0
                    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tadjusted VAF: " + str(vaf) + "\n")
        else:
            sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tselected VAF: " + str(vaf) + "\n")

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

        readtrack = dd(list)

        for readname in readlist:
            orig_name, readpos, pairend = readname.split(',')
            readtrack[orig_name].append('%s,%s' % (readpos, pairend))

        usedreads = 0
        newreadlist = []

        for orig_name in readtrack:
            for read_instance in readtrack[orig_name]:
                newreadlist.append(orig_name + ',' + read_instance)
                usedreads += 1

            if usedreads >= lastread:
                break

        readlist = newreadlist

        print "INFO\t" + now() + "\t" + mutid + "\tpicked: " + str(len(readlist)) + " reads"

        wrote = 0
        nmut = 0
        mut_out = {}
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
                mut_out[extqname] = read

        muts_written = {}

        for extqname in mut_out:
            if extqname not in muts_written:
                outbam_muts.write(mut_out[extqname])
                muts_written[extqname] = True

                if mutmates[extqname] is not None:
                    # is mate also in mutated list?
                    mate_read = mutmates[extqname]

                    pairname = 'F' # read is first in pair
                    if mate_read.is_read2:
                        pairname = 'S' # read is second in pair
                    if not mate_read.is_paired:
                        pairname = 'U' # read is unpaired

                    mateqname = ','.join((mate_read.qname,str(mate_read.pos),pairname))

                    if mateqname in mut_out:
                        # yes: output mutated mate
                        outbam_muts.write(mut_out[mateqname])
                        muts_written[mateqname] = True

                    else:
                        # no: output original mate
                        outbam_muts.write(mate_read)

        print "INFO\t" + now() + "\t" + mutid + "\twrote: " + str(wrote) + " reads, mutated: " + str(nmut) + " reads"

        if not hasSNP or args.force:
            outbam_muts.close()
            aligners.remap_bam(args.aligner, tmpoutbamname, args.refFasta, alignopts, mutid=mutid, paired=(not args.single), picardjar=args.picardjar)

            outbam_muts = pysam.Samfile(tmpoutbamname,'rb')
            coverwindow = 1
            incover  = countReadCoverage(bamfile,chrom,mutpos-coverwindow,mutpos+del_ln+coverwindow)
            outcover = countReadCoverage(outbam_muts,chrom,mutpos-coverwindow,mutpos+del_ln+coverwindow)

            avgincover  = float(sum(incover))/float(len(incover)) 
            avgoutcover = float(sum(outcover))/float(len(outcover))
            spikein_frac = 0.0
            if wrote > 0:
                spikein_frac = float(nmut)/float(wrote)

            # qc cutoff for final snv depth 
            if (avgoutcover > 0 and avgincover > 0 and avgoutcover/avgincover >= float(args.coverdiff)) or args.force:
                tmpbams.append(tmpoutbamname)
                indelstr = ''
                if is_insertion:
                    indelstr = ':'.join(('INS', chrom, str(start), ins))
                else:
                    indelstr = ':'.join(('DEL', chrom, str(start), str(end)))

                snvstr = chrom + ":" + str(start) + "-" + str(end) + " (VAF=" + str(vaf) + ")"
                log.write("\t".join(("indel",indelstr,str(mutpos),mutstr,str(avgincover),str(avgoutcover),str(spikein_frac),str(maxfrac)))+"\n")
            else:
                outbam_muts.close()
                os.remove(tmpoutbamname)
                if os.path.exists(tmpoutbamname + '.bai'):
                    os.remove(tmpoutbamname + '.bai')
                    
                print "WARN\t" + now() + "\t" + mutid + "\tdropped for outcover/incover < " + str(args.coverdiff)
                return None

        outbam_muts.close()
        bamfile.close()
        bammate.close()
        log.close() 

        return sorted(tmpbams)
        
    except Exception, e:
        sys.stderr.write("*"*60 + "\nencountered error in mutation spikein: " + mutid + "\n")
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

    aligners.checkoptions(args.aligner, alignopts, args.picardjar)

    # load readlist to avoid, if specified
    avoid = None
    if args.avoidreads is not None:
        avoid = dictlist(args.avoidreads)

    # make a temporary file to hold mutated reads
    outbam_mutsfile = "addindel." + str(uuid4()) + ".muts.bam"
    bamfile = pysam.Samfile(args.bamFileName, 'rb')
    outbam_muts = pysam.Samfile(outbam_mutsfile, 'wb', template=bamfile)
    outbam_muts.close()
    bamfile.close()
    tmpbams = []

    if not os.path.exists(args.tmpdir):
        os.mkdir(args.tmpdir)
        print "INFO\t" + now() + "\tcreated tmp directory: " + args.tmpdir

    if not os.path.exists('addindel_logs_' + os.path.basename(args.outBamFile)):
        os.mkdir('addindel_logs_' + os.path.basename(args.outBamFile))
        print "created directory: addindel_logs_" + os.path.basename(args.outBamFile)

    assert os.path.exists('addindel_logs_' + os.path.basename(args.outBamFile)), "could not create output directory!"
    assert os.path.exists(args.tmpdir), "could not create temporary directory!"

    pool = Pool(processes=int(args.procs))
    results = []

    ntried = 0
    for bedline in bedfile:
        if ntried < int(args.numsnvs) or int(args.numsnvs) == 0:
            c = bedline.strip().split()
            chrom = c[0]
            start = int(c[1])
            end   = int(c[2])
            vaf   = float(c[3])
            type  = c[4]
            ins   = None

            assert type in ('INS', 'DEL')
            if type == 'INS':
                ins = c[5]

            # make mutation (submit job to thread pool)
            result = pool.apply_async(makemut, [args, chrom, start, end, vaf, ins, avoid, alignopts])
            results.append(result)
            ntried += 1

    for result in results:
        try:
            tmpbamlist = result.get()
            if tmpbamlist is not None:
                for tmpbam in tmpbamlist:
                    if os.path.exists(tmpbam):
                        tmpbams.append(tmpbam)
        except AssertionError:
            print "****************************************************"
            print "* WARNING: assertion failed somewhere, check logs. *"
            print "****************************************************"

    if len(tmpbams) == 0:
        print "INFO\t" + now() + "\tno succesful mutations"
        sys.exit()

    tmpbams.sort()

    # merge tmp bams
    if len(tmpbams) == 1:
        os.rename(tmpbams[0],outbam_mutsfile)
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
        if args.tagreads:
            from bamsurgeon.markreads import markreads
            tmp_tag_bam = 'tag.%s.bam' % str(uuid4())
            markreads(outbam_mutsfile, tmp_tag_bam)
            move(tmp_tag_bam, outbam_mutsfile)
            print "INFO\t" + now() + "\ttagged reads."

        print "INFO\t" + now() + "\tdone making mutations, merging mutations into", args.bamFileName, "-->", args.outBamFile
        replace(args.bamFileName, outbam_mutsfile, args.outBamFile, seed=args.seed)

        #cleanup
        os.remove(outbam_mutsfile)
    
def run():
    # run this script
    parser = argparse.ArgumentParser(description='adds INDELs to reads, outputs modified reads as .bam along with mates')
    parser.add_argument('-v', '--varfile', dest='varFileName', required=True, help='Target regions to try and add a SNV, as BED')
    parser.add_argument('-f', '--bamfile', dest='bamFileName', required=True, help='sam/bam file from which to obtain reads')
    parser.add_argument('-r', '--reference', dest='refFasta', required=True, help='reference genome, fasta indexed with bwa index -a stdsw _and_ samtools faidx')
    parser.add_argument('-o', '--outbam', dest='outBamFile', required=True, help='.bam file name for output')
    parser.add_argument('-s', '--snvfrac', dest='snvfrac', default=1, help='maximum allowable linked SNP MAF (for avoiding haplotypes) (default = 1)')
    parser.add_argument('-m', '--mutfrac', dest='mutfrac', default=0.5, help='allelic fraction at which to make SNVs (default = 0.5)')
    parser.add_argument('-n', '--numsnvs', dest='numsnvs', default=0, help="maximum number of mutations to try (default: entire input)")
    parser.add_argument('-c', '--cnvfile', dest='cnvfile', default=None, help="tabix-indexed list of genome-wide absolute copy number values (e.g. 2 alleles = no change)")
    parser.add_argument('-d', '--coverdiff', dest='coverdiff', default=0.1, help="allow difference in input and output coverage (default=0.1)")
    parser.add_argument('-p', '--procs', dest='procs', default=1, help="split into multiple processes (default=1)")
    parser.add_argument('--picardjar', default=None, help='path to picard.jar')
    parser.add_argument('--mindepth', default=10, help='minimum read depth to make mutation (default = 10)')
    parser.add_argument('--maxdepth', default=2000, help='maximum read depth to make mutation (default = 2000)')
    parser.add_argument('--minmutreads', default=3, help='minimum number of mutated reads to output per site')
    parser.add_argument('--avoidreads', default=None, help='file of read names to avoid (mutations will be skipped if overlap)')
    parser.add_argument('--nomut', action='store_true', default=False, help="dry run")
    parser.add_argument('--det', action='store_true', default=False, help="deterministic base changes: make transitions only")
    parser.add_argument('--force', action='store_true', default=False, help="force mutation to happen regardless of nearby SNP or low coverage")
    parser.add_argument('--single', action='store_true', default=False, help="input BAM is simgle-ended (default is paired-end)")
    parser.add_argument('--maxopen', dest='maxopen', default=1000, help="maximum number of open files during merge (default 1000)")
    parser.add_argument('--requirepaired', action='store_true', default=False, help='skip mutations if unpaired reads are present')
    parser.add_argument('--aligner', default='backtrack', help='supported aligners: ' + ','.join(aligners.supported_aligners_bam))
    parser.add_argument('--alignopts', default=None, help='aligner-specific options as comma delimited list of option1:value1,option2:value2,...')
    parser.add_argument('--tagreads', action='store_true', default=False, help='add BS tag to altered reads')
    parser.add_argument('--skipmerge', action='store_true', default=False, help="final output is tmp file to be merged")
    parser.add_argument('--ignorepileup', action='store_true', default=False, help="do not check pileup depth in mutation regions")
    parser.add_argument('--tmpdir', default='addindel.tmp', help='temporary directory (default=addindel.tmp)')
    parser.add_argument('--seed', default=None, help='seed random number generation')
    args = parser.parse_args()
    main(args)

if __name__ == '__main__':
    run()
