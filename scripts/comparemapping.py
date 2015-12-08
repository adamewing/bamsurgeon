#!/home/tmlundbe/anaconda2/bin/python

import subprocess,pysam, sys, re

def namesortbam(inbam,outbam):
    print "sorting by name:",inbam,"-->",outbam
    outbam = re.sub('.bam$','',outbam)

    sortargs = ['samtools','sort','-n','-m','20000000000',inbam, '-o', outbam]
    subprocess.call(sortargs)

def compare(bamfile1,bamfile2):
    bam1sort = re.sub('.bam$','',bamfile1) + ".namesort.bam"
    bam2sort = re.sub('.bam$','',bamfile2) + ".namesort.bam"

    namesortbam(sys.argv[1],bam1sort)
    namesortbam(sys.argv[2],bam2sort)

    bam1 = pysam.Samfile(bam1sort,'rb')
    bam2 = pysam.Samfile(bam2sort,'rb')

    n = {}
    n['total'] = 0 # read count
    n['diff']  = 0 # number of reads where mapping coords change
    n['unmap'] = 0 # number of reads where mappable-->unmappable
    n['map']   = 0 # number of reads where unmappable-->mappable
    n['clip']  = 0 # number of reads where clipping changes
    n['rep']   = 0 # number of read pairs with a repeat (mapping qual 0) read
    n['seq']   = 0 # number of reads with a sequence change 

    for read1 in bam1.fetch(until_eof=True):
        read2 = bam2.next()
        assert read1.qname == read2.qname # require identical sets of read names

        if read1.seq != read2.seq:
            n['seq'] += 1

        if read1.mapq > 1 and read2.mapq > 1:
            n['total'] += 1
            if not read1.is_unmapped and read2.is_unmapped:
                n['unmap'] += 1
            elif read1.is_unmapped and not read2.is_unmapped:
                n['map'] += 1
            elif not read1.is_unmapped and not read2.is_unmapped:
                if read1.alen != read2.alen:
                    n['clip'] += 1
                if read1.pos != read2.pos:
                    n['diff'] += 1
        else:
            n['rep'] += 1

    bam1.close()
    bam2.close()

    for stat,val in n.iteritems():
        print stat,val

if len(sys.argv) != 3:
    print "this script compares two .bams with identical sets of read names (required) and reports mapping differences"
    print "usage:",sys.argv[0],"<bamfile1.bam> <bamfile2.bam>"
else:
    compare(sys.argv[1],sys.argv[2])
