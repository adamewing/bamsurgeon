#!/usr/bin/env python

import argparse
import os
import subprocess


def getsegs(bam, mindepth, minlength):
    seglist = []
    cmd = ['samtools', 'mpileup', bam]
    
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)

    seg = {'chrom': None, 'start': None, 'end': None}

    for line in p.stdout:
        chrom, pos, base, depth = line.strip().split()[:4]

        depth = int(depth)
        pos   = int(pos)

        if seg['chrom'] != chrom:
            if seg['chrom'] is not None:
                if seg['end'] - seg['start'] >= minlength:
                    seglist.append(seg)

            seg = {'chrom': chrom, 'start': pos, 'end': pos}

        else:
            if pos == seg['end']+1 and depth > mindepth:
                seg['end'] = pos

            else:
                if seg['end'] - seg['start'] >= minlength:
                    seglist.append(seg)
                seg = {'chrom': chrom, 'start': pos, 'end': pos}

    return seglist


def main(args):
    assert os.path.exists(args.bam + '.bai'), ".bai index not found for BAM: " + args.bam

    seglist = getsegs(args.bam, int(args.depth), int(args.length))
    for seg in seglist:
        print('\t'.join((seg['chrom'], str(seg['start']), str(seg['end']))))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="return segments of the genome covered at some minimum depth")
    parser.add_argument('-b', '--bam', required=True, help="indexed BAM file")
    parser.add_argument('-d', '--depth', default=10, help="minimum depth to report segment (default=10)")
    parser.add_argument('-l', '--length', default=1000, help="minimum length covered segment to report (default=1000)")

    args = parser.parse_args()
    main(args)
