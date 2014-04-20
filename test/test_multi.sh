#!/bin/sh

# adds up to 100 SNPs to a ~770 kb region around the LARGE gene
# requires samtools/bcftools

if [ $# -ne 1 ]
then
    echo "usage: $0 <reference indexed with bwa index>"
    exit 65
fi

if [ ! -e ../addsnv.py ]
then
    echo "addsnv.py isn't one directory level down (../addsnv.py) as expected"
    exit 65
fi

if [ ! -e $1 ]
then
    echo "can't find reference .fasta: $1, please supply a bwa-indexed .fasta"
    exit 65
fi

if [ ! -e $1.bwt ]
then
    echo "can't find $1.bwt: is $1 indexed with bwa?"
    exit 65
fi

../addsv.py -p 4 -v ../test_data/test_multibam_sv.txt -f ../test_data/test_multi.lib1.bam:../test_data/test_multi.lib2.bam:../test_data/test_multi.lib3.bam -r $1 -o ../test_data/test_multi.lib1.sv.bam:../test_data/test_multi.lib2.sv.bam:../test_data/test_multi.lib3.sv.bam -c ../test_data/test_multi.cnv.txt.gz
if [ $? -ne 0 ]
then
  echo "addsv.py failed."
  exit 65
else
  echo "sorting output bam..."
  samtools sort ../test_data/test_multi.lib1.sv.bam ../test_data/test_multi.lib1.sv.bam.sorted
  samtools sort ../test_data/test_multi.lib2.sv.bam ../test_data/test_multi.lib2.sv.bam.sorted
  samtools sort ../test_data/test_multi.lib3.sv.bam ../test_data/test_multi.lib3.sv.bam.sorted
  
  mv ../test_data/test_multi.lib1.sv.bam.sorted ../test_data/test_multi.lib1.sv.bam 
  mv ../test_data/test_multi.lib2.sv.bam.sorted ../test_data/test_multi.lib2.sv.bam 
  mv ../test_data/test_multi.lib3.sv.bam.sorted ../test_data/test_multi.lib3.sv.bam 

  echo "indexing output bam..."
  samtools index ../test_data/test_multi.lib1.sv.bam
  samtools index ../test_data/test_multi.lib2.sv.bam
  samtools index ../test_data/test_multi.lib3.sv.bam
  
  echo "making pileups..."
  samtools mpileup $2 ../test_data/test_multi.lib1.bam ../test_data/test_multi.lib2.bam ../test_data/test_multi.lib3.bam ../test_data/test_multi.lib1.sv.bam ../test_data/test_multi.lib2.sv.bam ../test_data/test_multi.lib3.sv.bam > test_multibam_sv.pileup.txt
  echo "done. output in test_sv.pileup.txt"
fi
