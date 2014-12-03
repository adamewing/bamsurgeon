#!/bin/bash

# adds up to 100 SNPs to a ~770 kb region around the LARGE gene
# requires samtools/bcftools

if [ $# -ne 6 ]
then
    echo "usage: $0 <number of SNPs> <number of threads> <reference indexed with samtools> <path to SamToFastq.jar provided by picard tools> <gsnap reference dir> <gsnap reference name>"
    exit 65
fi

if [ ! -e ../addsnv.py ]
then
    echo "addsnv.py isn't one directory level down (../addsnv.py) as expected"
    exit 65
fi

if [ ! -e $3 ]
then
    echo "can't find reference .fasta: $3, please supply a samtools-indexed .fasta"
    exit 65
fi

if [ ! -e $4 ]
then
    echo "cannot find SamToFastq.jar"
    exit 65
fi

if [ ! -e $5 ]
then
    echo "cannot find gsnap ref dir"
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

../addsnv.py -v ../test_data/random_snvs.txt -f ../test_data/testregion_gsnap.bam -r $3 -o ../test_data/testregion_gsnap_mut.bam -n $1 -c ../test_data/test_cnvlist.txt.gz -p $2 --ignoresnps --samtofastq $4 --aligner gsnap --alignopts gsnaprefdir:$5,gsnaprefname:$6
if [ $? -ne 0 ]
then
 echo "addsnv.py failed. Are all the prequisites installed?"
 exit 65
else
    samtools sort ../test_data/testregion_gsnap_mut.bam ../test_data/testregion_gsnap_mut.sorted
    mv ../test_data/testregion_gsnap_mut.sorted.bam ../test_data/testregion_gsnap_mut.bam
    samtools index ../test_data/testregion_gsnap_mut.bam
    samtools mpileup -ugf $3 ../test_data/testregion_gsnap_mut.bam | bcftools view -bvcg - > result.raw.bcf
    bcftools view result.raw.bcf
fi
