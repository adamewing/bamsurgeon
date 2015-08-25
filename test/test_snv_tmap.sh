#!/bin/bash

# adds up to 100 SNPs to a ~770 kb region around the LARGE gene
# requires samtools/bcftools

if [ $# -ne 4 ]
then
    echo "usage: $0 <number of SNPs> <number of threads> <reference indexed with tmap index> <path to picard.jar>"
    exit 65
fi

if [ ! -e ../addsnv.py ]
then
    echo "addsnv.py isn't one directory level down (../addsnv.py) as expected"
    exit 65
fi

if [ ! -e $4 ]
then
    echo "cannot find SamToFastq.jar"
    exit 65
fi

if [ ! -e $3 ]
then
    echo "can't find reference .fasta: $3, please supply a bwa-indexed .fasta"
    exit 65
fi

if [ ! -e $3.bwt ]
then
    echo "can't find $3.bwt: is $3 indexed with bwa?"
    exit 65
fi

if ! [[ $1 =~ ^[0-9]+$ ]]
then
    echo "arg 1 must be an integer (number of SNVs to add)"
    exit 65
fi

if [ $1 -gt 100 ]
then
    echo "max number of SNVs must be <= 100"
    exit 65
fi

if [ $1 -lt 1 ]
then
    echo "min number of SNVs must be > 0"
    exit 65
fi

../addsnv.py --single -v ../test_data/random_snvs.txt -f ../test_data/testregion_tmap.bam -r $3 -o ../test_data/testregion_mut.bam -n $1 -p $2 --picardjar $4 --aligner tmap 

if [ $? -ne 0 ]
then
 echo "addsnv.py failed. Are all the prequisites installed?"
 exit 65
else
  samtools mpileup -ugf $3 ../test_data/testregion_mut.bam | bcftools call -vm
fi
