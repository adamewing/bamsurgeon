#!/usr/bin/env python


import subprocess
import pysam
import random
import hashlib
import string
import os

from collections import Counter
from shutil import move
from re import sub

import logging
FORMAT = '%(levelname)s %(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

random_salt = None


def read_hash_fraction(query_name):
    global random_salt
    if random_salt is None:
        random_salt = ''.join(random.choice(string.ascii_lowercase + string.ascii_uppercase + string.digits) for i in range(10))
    read_hash = int(hashlib.md5((random_salt + query_name).encode()).hexdigest(), 16)
    read_random_factor = (read_hash % 1000000) / 1000000.0
    return read_random_factor


def get_avg_coverage(alignment_file, chrom, start, end):
    split_coverage = alignment_file.count_coverage(chrom, start, end, quality_threshold=0)
    base_sum = [sum(x) for x in split_coverage]
    return sum(base_sum) / float(end - start)


def rc(dna):
    ''' reverse complement '''
    complements = str.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
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
    logger.info("len(bamlist): %d" % len(bamlist))

    if len(bamlist) == 1:
        logger.info("only one BAM to merge, renaming " + bamlist[0] + " --> " + outbamfn)
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

        for submergefn, tmpbams in merge_sublists.items():
            if len(tmpbams) == 1:
                move(tmpbams[0], submergefn)
                logger.info("renamed: " + tmpbams[0] + " --> " + submergefn)
            else:
                args = ['samtools', 'merge', '-f', submergefn] + tmpbams
                logger.info("merging, cmd: " + ' '.join(args))
                subprocess.check_call(args)

        if len(merge_sublists.keys()) == 1:
            logger.info("merge finished, renaming: " + list(merge_sublists.keys())[0] + " --> " + outbamfn)
            move(list(merge_sublists.keys())[0], outbamfn)
        else:
            args = ['samtools', 'merge', '-f', outbamfn] + list(merge_sublists.keys())
            logger.info("final merge, cmd: " + ' '.join(args))
            subprocess.check_call(args)

        for submergefn in merge_sublists.keys():
            if os.path.exists(submergefn):
                os.remove(submergefn)

    if not debug:
        for bamfile in bamlist:
            if os.path.exists(bamfile):
                os.remove(bamfile)
            if os.path.exists(bamfile + '.bai'):
                os.remove(bamfile + '.bai')


def double_fastq_to_interleaved(fq1, fq2, outfq):
    def get_fastq_read_iterator(fq):
        f = open(fq, 'r')
        while True:
            name = f.readline()
            if not name:
                break
            other = [f.readline() for i in range(3)]
            yield [name] + other
        f.close()
    f_1_iter = get_fastq_read_iterator(fq1)
    f_2_iter = get_fastq_read_iterator(fq2)
    f_out = open(outfq, 'w')
    while True:
        try:
            f_out.write(''.join(next(f_1_iter)))
            f_out.write(''.join(next(f_2_iter)))
        except StopIteration:
            break
    f_out.close()


def bamtofastq(bam, picardjar, threads=1, paired=True, twofastq=False):
    ''' if twofastq is True output two fastq files instead of interleaved (default) for paired-end'''
    assert bam.endswith('.bam')

    outfq = None
    outfq_pair = None

    cmd = ['java', '-XX:ParallelGCThreads=' + str(threads), '-jar', picardjar, 'SamToFastq', 'VALIDATION_STRINGENCY=SILENT', 'INPUT=' + bam]
    cmd.append('INCLUDE_NON_PRIMARY_ALIGNMENTS=false')  # in case the default ever changes

    logger.info("converting BAM " + bam + " to FASTQ\n")
    if paired:
        if twofastq:  # two-fastq paired end
            outfq_pair = [sub('bam$', '1.fastq', bam), sub('bam$', '2.fastq', bam)]
            cmd.append('F=' + outfq_pair[0])
            cmd.append('F2=' + outfq_pair[1])
            subprocess.check_call(cmd)
            assert os.path.exists(outfq_pair[0]) and os.path.exists(outfq_pair[1])
            return outfq_pair
        else:  # interleaved paired-end
            outfq = sub('bam$', 'fastq', bam)
            cmd.append('FASTQ=' + outfq)
            cmd.append('INTERLEAVE=true')
            subprocess.check_call(cmd)
            assert os.path.exists(outfq)  # conversion failed
            return [outfq]
    else:
        outfq = sub('bam$', 'fastq', bam)
        cmd.append('FASTQ=' + outfq)
        subprocess.check_call(cmd)
        return [outfq]


def check_min_read_count(bamfile, fasta_ref, min_num_reads=0):
    # Returns True if the BAM has at least min_num_reads reads
    bam = pysam.AlignmentFile(bamfile, reference_filename=fasta_ref)
    # If the BAM is indexed, we can just check the header
    try:
        if bam.is_bam and bam.check_index():
            return bam.mapped + bam.unmapped > min_num_reads
    except AttributeError:
        pass
    except ValueError:
        pass
    # If the BAM is not indexed, we have to count the reads manually
    for i, _ in enumerate(bam.fetch(until_eof=True)):
        if i + 1 > min_num_reads:
            return True
    return False


def dictlist(fn):
    d = {}
    with open(fn, 'r') as inlist:
        for name in inlist:
            d[name.strip()] = True
    return d
