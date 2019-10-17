#!/bin/sh

# requires samtools/bcftools

if [ $# -ne 1 ]
then
    echo "usage: $0 <path to picard.jar>"
    exit 65
fi

REF=../test_data/Homo_sapiens_chr22_assembly19.fasta

command -v addindel.py >/dev/null 2>&1 || { echo "addindel.py isn't installed" >&2; exit 65; }

if [ ! -e $1 ]
then
    echo "cannot find picard.jar"
    exit 65
fi

export BAMSURGEON_PICARD_JAR=$1

addindel.py -v ../test_data/test_indels.txt -f  ../test_data/testregion_realign.bam -r $REF -o ../test_data/testregion_mut.bam  --aligner mem --seed 1234

if [ $? -ne 0 ]
then
 echo "addindel.py failed."
 exit 1 
else
    samtools sort -T ../test_data/testregion_mut.sorted.bam -o ../test_data/testregion_mut.sorted.bam ../test_data/testregion_mut.bam
    mv ../test_data/testregion_mut.sorted.bam ../test_data/testregion_mut.bam
    samtools index ../test_data/testregion_mut.bam
    samtools mpileup -ugf $REF ../test_data/testregion_mut.bam | bcftools call -vm
fi
