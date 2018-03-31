#!/usr/bin/env python


import datetime
import subprocess
import pysam
import os
import sys

from string import maketrans
from collections import Counter
from shutil import move
from re import sub


def now():
    return str(datetime.datetime.now())


def rc(dna):
    ''' reverse complement '''
    complements = maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
    return dna.translate(complements)[::-1]


def majorbase(basepile):
    """returns tuple: (major base, count)
    """
    return Counter(basepile).most_common()[0]


def minorbase(basepile):
    """returns tuple: (minor base, count)
    """
    c = Counter(basepile)
    if len(list(c.elements())) > 1:
        return c.most_common(2)[-1]
    else:
        return c.most_common()[0]


def mergebams(bamlist, outbamfn, maxopen=100, debug=False):
    ''' call samtools to merge a list of bams hierarchically '''

    assert outbamfn.endswith('.bam')
    print "INFO\t" + now() + "\tlen(bamlist)", len(bamlist)

    if len(bamlist) == 1:
        print "INFO\t" + now() + "\tonly one BAM to merge, renaming",bamlist[0],"-->",outbamfn
        move(bamlist[0], outbamfn)
    else:
        nmerge = 1
        mergenum = 0
        merge_sublists = {}
        for tmpbam in bamlist:
            mergetmp = "tmp.merging." + str(mergenum) + "." + outbamfn
            if mergetmp in merge_sublists:
                merge_sublists[mergetmp].append(tmpbam)
            else:
                merge_sublists[mergetmp] = []
                merge_sublists[mergetmp].append(tmpbam)
            if nmerge % maxopen == 0:
                mergenum += 1
            nmerge += 1

        for submergefn, tmpbams in merge_sublists.iteritems():
            if len(tmpbams) == 1:
                move(tmpbams[0], submergefn)
                print "INFO\t" + now() + "\trenamed:",tmpbams[0], " --> ", submergefn
            else:
                args = ['samtools','merge','-f',submergefn] + tmpbams 
                print "INFO\t" + now() + "\tmerging, cmd: ", args
                subprocess.call(args)

        if len(merge_sublists.keys()) == 1:
            print "INFO\t" + now() + "\tmerge finished, renaming: ", merge_sublists.keys()[0]," --> ", outbamfn
            move(merge_sublists.keys()[0], outbamfn)
        else:
            args = ['samtools','merge','-f',outbamfn] + merge_sublists.keys()
            print "INFO\t" + now() + "\tfinal merge, cmd: ", args
            subprocess.call(args)

        for submergefn in merge_sublists.keys():
            if os.path.exists(submergefn):
                os.remove(submergefn)

    if not debug:
        for bamfile in bamlist:
            if os.path.exists(bamfile):
                os.remove(bamfile)
            if os.path.exists(bamfile + '.bai'):
                os.remove(bamfile + '.bai')


def bamtofastq(bam, picardjar, threads=1, paired=True, twofastq=False):
    ''' if twofastq is True output two fastq files instead of interleaved (default) for paired-end'''
    assert os.path.exists(picardjar)
    assert bam.endswith('.bam')

    outfq = None
    outfq_pair = None

    cmd = ['java', '-XX:ParallelGCThreads=' + str(threads), '-jar', picardjar, 'SamToFastq', 'VALIDATION_STRINGENCY=SILENT', 'INPUT=' + bam]
    cmd.append('INCLUDE_NON_PRIMARY_ALIGNMENTS=false') # in case the default ever changes
    
    if paired:
        if twofastq: # two-fastq paired end
            outfq_pair = [sub('bam$', '1.fastq', bam), sub('bam$', '2.fastq', bam)]
            cmd.append('F=' + outfq_pair[0])
            cmd.append('F2=' + outfq_pair[1])
        else: # interleaved paired-end
            outfq = sub('bam$', 'fastq', bam)
            cmd.append('FASTQ=' + outfq)
            cmd.append('INTERLEAVE=true')
    else:
        outfq = sub('bam$', 'fastq', bam)
        cmd.append('FASTQ=' + outfq)

    sys.stdout.write("INFO\t" + now() + "\tconverting BAM " + bam + " to FASTQ\n")
    subprocess.call(cmd)

    if outfq is not None:
        assert os.path.exists(outfq) # conversion failed
        return [outfq]

    if outfq_pair is not None:
        assert os.path.exists(outfq_pair[0]) and os.path.exists(outfq_pair[1])
        return outfq_pair

    return None

    
def fastqreadcount(fastqfile):
    assert not fastqfile.endswith('gz') # not supported yet
    return sum(1 for line in open(fastqfile))/4


def bamreadcount(bamfile):
    bam = pysam.Samfile(bamfile, 'rb')
    if os.path.exists(bamfile + '.bai'):
        return bam.mapped + bam.unmapped
    else:
        return(list(bam.fetch(until_eof=True)))


def dictlist(fn):
    d = {}
    with open(fn, 'r') as inlist:
        for name in inlist:
            d[name.strip()] = True
    return d
