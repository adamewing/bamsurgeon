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


addindel.py -v ../test_data/test_ont_indel_ins.txt -f ../test_data/testregion_ont.bam -r $REF -o ../test_data/testregion_ont_mut.bam -n 5 --picardjar $1 --aligner minimap2 --alignopts x:map-ont --seed 1234 --ignorepileup --single --mindepth 4

samtools sort ../test_data/testregion_ont_mut.bam > ../test_data/testregion_ont_mut.sorted.bam
mv ../test_data/testregion_ont_mut.sorted.bam ../test_data/testregion_ont_mut.bam
samtools index ../test_data/testregion_ont_mut.bam

#if [ $? -ne 0 ]
#then
# echo "addsnv.py failed. Are all the prequisites installed?"
# exit 1 
#else
#  samtools mpileup -ugf $REF ../test_data/testregion_mut.bam | bcftools call -vm
#fi
