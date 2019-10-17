#!/bin/bash

if [ $# -ne 1 ]
then
    echo "usage: $0 <path to picard.jar>"
    exit 65
fi

REF=../test_data/Homo_sapiens_chr22_assembly19.fasta

command -v addsnv.py >/dev/null 2>&1 || { echo "addsnv.py isn't installed" >&2; exit 65; }

if [ ! -e $1 ]
then
    echo "cannot find picard.jar"
    exit 65
fi

export BAMSURGEON_PICARD_JAR=$1

addsnv.py -v ../test_data/random_snvs.txt -f ../test_data/testregion_realign.bam -r $REF -o ../test_data/testregion_mut.bam -n 5 -c ../test_data/test_cnvlist.txt.gz  --aligner mem --seed 1234

if [ $? -ne 0 ]
then
 echo "addsnv.py failed. Are all the prequisites installed?"
 exit 1 
else
  samtools mpileup -ugf $REF ../test_data/testregion_mut.bam | bcftools call -vm
fi
