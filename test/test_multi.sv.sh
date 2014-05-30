#!/bin/sh

# adds up to 100 SNPs to a ~770 kb region around the LARGE gene
# requires samtools/bcftools

if [ $# -ne 1 ]
then
    echo "usage: $0 <reference indexed with bwa index>"
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

addsv -p 2 -v ../test_data/test_multi.sv.txt -t ../test_data/test_multi.sv.pct -i 1_259_22:1_371_28:0.4_540_50,0.6_143_50 -f ../test_data/test_multi.lib1.bam:../test_data/test_multi.lib2.bam:../test_data/test_multi.lib3.bam -r $1 -o ../test_data/test_multi.sv.lib1.bam:../test_data/test_multi.sv.lib2.bam:../test_data/test_multi.sv.lib3.bam
if [ $? -ne 0 ]
then
  echo "addsv.py failed."
  exit 65
else
  echo "sorting output bam..."
  samtools sort ../test_data/test_multi.sv.lib1.bam ../test_data/test_multi.sv.lib1.bam.sorted
  samtools sort ../test_data/test_multi.sv.lib2.bam ../test_data/test_multi.sv.lib2.bam.sorted
  samtools sort ../test_data/test_multi.sv.lib3.bam ../test_data/test_multi.sv.lib3.bam.sorted
  
  mv ../test_data/test_multi.sv.lib1.bam.sorted.bam ../test_data/test_multi.sv.lib1.bam 
  mv ../test_data/test_multi.sv.lib2.bam.sorted.bam ../test_data/test_multi.sv.lib2.bam 
  mv ../test_data/test_multi.sv.lib3.bam.sorted.bam ../test_data/test_multi.sv.lib3.bam 

  echo "indexing output bam..."
  samtools index ../test_data/test_multi.sv.lib1.bam
  samtools index ../test_data/test_multi.sv.lib2.bam
  samtools index ../test_data/test_multi.sv.lib3.bam
  
  echo "making pileups..."
  samtools mpileup $2 ../test_data/test_multi.lib1.bam ../test_data/test_multi.lib2.bam ../test_data/test_multi.lib3.bam ../test_data/test_multi.sv.lib1.bam ../test_data/test_multi.sv.lib2.bam ../test_data/test_multi.sv.lib3.bam > ../test_data/test_multi.svbam_sv.pileup.txt
  echo "done. output in test_sv.pileup.txt"
  
  echo "collecting mutations..."
  cat addsv_logs_test_multi.sv.lib1.bam/*.log | grep -E "dup|del|ins|inv" >test_multi.sv.tab

fi
