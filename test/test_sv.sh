#!/bin/bash

# adds up to 100 SNPs to a ~770 kb region around the LARGE gene
# requires samtools/bcftools

if [ $# -ne 2 ]
then
    echo "usage: $0 <number of threads> <reference indexed with bwa index>"
    exit 65
fi

command -v addsv.py >/dev/null 2>&1 || { echo "addsv.py isn't installed" >&2; exit 65; }

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

if [ ! -e $2 ]
then
    echo "can't find reference .fasta: $2, please supply a bwa-indexed .fasta"
    exit 65
fi

if [ ! -e $2.bwt ]
then
    echo "can't find $2.bwt: is $2 indexed with bwa?"
    exit 65
fi

addsv.py -p $1 -v ../test_data/test_sv.txt -f ../test_data/testregion.bam -r $2 -o ../test_data/testregion_sv_mut.bam -c ../test_data/test_cnvlist.txt.gz --seed 1234
if [ $? -ne 0 ]
then
  echo "addsv.py failed."
  exit 65
else
  echo "sorting output bam..."
  samtools sort ../test_data/testregion_sv_mut.bam ../test_data/testregion_sv_mut.sorted
  mv ../test_data/testregion_sv_mut.sorted.bam ../test_data/testregion_sv_mut.bam

  echo "indexing output bam..."
  samtools index ../test_data/testregion_sv_mut.bam

  echo "making pileups..."
  samtools mpileup -f $2 ../test_data/testregion_sv_mut.bam ../test_data/testregion.bam > test_sv.pileup.txt
  echo "done. output in test_sv.pileup.txt"
fi
