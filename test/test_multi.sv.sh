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

addsv.py -p 2 -v spikein.txt -i 1_259_22:1_371_28:0.4_540_50,0.6_143_50 -f VU000123.chr11.65M-75M.lib1.bam:VU000123.chr11.65M-75M.lib2.bam:VU000123.chr11.65M-75M.lib3.bam -r $1 -o VU000123.chr11.65M-75M.lib1.sv.bam:VU000123.chr11.65M-75M.lib2.sv.bam:VU000123.chr11.65M-75M.lib3.sv.bam
if [ $? -ne 0 ]
then
  echo "addsv.py failed."
  exit 65
else
  echo "sorting output bam..."
  samtools sort VU000123.chr11.65M-75M.lib1.sv.bam VU000123.chr11.65M-75M.lib1.sv.bam.sorted
  samtools sort VU000123.chr11.65M-75M.lib2.sv.bam VU000123.chr11.65M-75M.lib2.sv.bam.sorted
  samtools sort VU000123.chr11.65M-75M.lib3.sv.bam VU000123.chr11.65M-75M.lib3.sv.bam.sorted
  
  mv VU000123.chr11.65M-75M.lib1.sv.bam.sorted.bam VU000123.chr11.65M-75M.lib1.sv.bam 
  mv VU000123.chr11.65M-75M.lib2.sv.bam.sorted.bam VU000123.chr11.65M-75M.lib2.sv.bam 
  mv VU000123.chr11.65M-75M.lib3.sv.bam.sorted.bam VU000123.chr11.65M-75M.lib3.sv.bam 

  echo "indexing output bam..."
  samtools index VU000123.chr11.65M-75M.lib1.sv.bam
  samtools index VU000123.chr11.65M-75M.lib2.sv.bam
  samtools index VU000123.chr11.65M-75M.lib3.sv.bam
  
  echo "making pileups..."
  samtools mpileup $2 VU000123.chr11.65M-75M.lib1.bam VU000123.chr11.65M-75M.lib2.bam VU000123.chr11.65M-75M.lib3.bam VU000123.chr11.65M-75M.lib1.sv.bam VU000123.chr11.65M-75M.lib2.sv.bam VU000123.chr11.65M-75M.lib3.sv.bam > VU000123.chr11.65M-75Mbam_sv.pileup.txt
  echo "done. output in test_sv.pileup.txt"
fi
