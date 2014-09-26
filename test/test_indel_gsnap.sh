#!/bin/sh

# requires samtools/bcftools

if [ $# -ne 5 ]
then
    echo "usage: $0 <number of threads> <reference indexed with bwa index> <path to SamToFastq.jar provided by picard tools> <gsnap reference dir> <gsnap reference name>"
    exit 65
fi

if [ ! -e ../addindel.py ]
then
    echo "addsnv.py isn't one directory level down (../addsnv.py) as expected"
    exit 65
fi

if [ ! -e $3 ]
then
    echo "cannot find SamToFastq.jar"
    exit 65
fi

if [ ! -e $4 ]
then
    echo "cannot find gsnap ref dir"
    exit 65
fi

if [ ! -e $2 ]
then
    echo "can't find reference .fasta: $2, please supply a samtools-indexed .fasta"
    exit 65
fi

if [ ! -e $2.fai ]
then
    echo "can't find $2.fai: is $2 indexed with samtools faidx?"
    exit 65
fi


../addindel.py -v ../test_data/test_indels.txt -f ../test_data/testregion_gsnap.bam -r $2 -o ../test_data/testregion_gsnap_mut.bam -p $1 --samtofastq $3 --gsnap --gsnaprefdir $4 --gsnaprefname $5
if [ $? -ne 0 ]
then
 echo "addindel.py failed."
 exit 65
else
    samtools sort ../test_data/testregion_gsnap_mut.bam ../test_data/testregion_gsnap_mut.sorted
    mv ../test_data/testregion_gsnap_mut.sorted.bam ../test_data/testregion_gsnap_mut.bam
    samtools index ../test_data/testregion_gsnap_mut.bam
    samtools mpileup -ugf $2 ../test_data/testregion_gsnap_mut.bam | bcftools view -bvcg - > result.raw.bcf
    bcftools view result.raw.bcf
fi
