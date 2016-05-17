#!/usr/bin/env python

import pysam
import sys
from re import sub
from random import random
import math

# GLOBALs
global readcount
global outbams_reads        # hash of outbam -> number reads written to it, every window
global bam_percents         # hash of outbam -> percents
global total_reads          # hash of outbam -> total number of reads written to file


# HELPERs
def error_check (outbams, percents):
    percentSum = sum(float(i) for i in percents)

    assert len(outbams) == len(percents), "List of output BAMs not same length as list of percents"
    assert abs(percentSum-1) < 1e-05 , "Percents don't add to 1 +- 1e-05"

    # safety measures that account for remaining < 1e-05
    if percentSum != 1:
        percents[len(percents)-1] = 1 - sum(float(percents[i]) for i in range(0,len(percents)-1))
    return percents

def reset ():
    global readcount
    global outbams_reads

    readcount = 0
    for outbam in outbams_reads:
        try:
            total_reads[outbam] += outbams_reads[outbam]
        except:
            total_reads[outbam] = 0
        outbams_reads[outbam] = 0


def findBam (rnd):
    # remaining outbams in current window to distribute remaining reads
    outbams_left = []
    # aggregated bamsplit percentages
    percents_left = []

    for outbam in outbams_reads:
        if outbams_reads[outbam] != int(math.ceil(float(bam_percents[outbam]) * window)):
            try:
                val = float(percents_left[len(percents_left)-1]) + float(bam_percents[outbam])
            except:
                val = float(bam_percents[outbam])
            outbams_left.append(outbam)
            percents_left.append(val)

    assert(len(percents_left) == len(outbams_left))

    idx = 0
    while True:
        if rnd <= float(percents_left[idx]/max(percents_left)) or idx == len(percents_left)-1:
            outbam = outbams_left[idx]
            outbams_reads[outbam] += 1
            return outbam
        idx += 1
    
    print "Invalid random value: ", rnd
    sys.exit(1)


# MAIN
# initialize
readcount = 0
outbams_reads = {}
total_reads = {}
bam_percents = {}
window = 1000


if len(sys.argv) == 5:
    assert sys.argv[1].endswith('.bam')

    outdir = sys.argv[2]
    outbams = sys.argv[3].split(',')
    percents = sys.argv[4].split(',')

    percents = error_check(outbams, percents)

    # file handles
    inbam_fh = pysam.Samfile(sys.argv[1], 'rb')
    outbams_fh = {}

    all_reads = {}              # hash of read name -> outbam

    # initialize some hashes
    for idx in range(len(outbams)):
        bam_percents[outbams[idx]] = percents[idx]
        outbams_fh[outbams[idx]] = pysam.Samfile(outdir+'/'+outbams[idx], 'wb', template=inbam_fh)
        outbams_reads[outbams[idx]] = 0

    reset()
   
    for read in inbam_fh.fetch(until_eof=True):
        if not read.is_secondary:
            # if its read already seen
            if read.qname in all_reads:
                outbam = all_reads[read.qname]

            else:
                readcount += 1
                rnd = random()
            
                outbam = findBam(rnd)
                all_reads[read.qname] = outbam

            outbam_fh = outbams_fh[outbam]
            outbam_fh.write(read)

            if readcount == window:
                reset()

        elif read.is_secondary:
            print "Skipping read: ", read.qname

    inbam_fh.close()
    for outbam in outbams_fh:
        outbams_fh[outbam].close()

    reset()
    print "TOTAL READS: ", total_reads
else:
    print "usage:",sys.argv[0],"<in.bam> <output_dir> <list of out.bams> <list of percentage splits>"
