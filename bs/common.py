#!/usr/bin/env python


import datetime
import subprocess
import os

from shutil import move


def now():
    return str(datetime.datetime.now())


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
                os.remove(bamfile + '.bai')


