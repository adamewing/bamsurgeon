#!/bin/sh

# adds up to 100 SNPs to a ~770 kb region around the LARGE gene
# requires samtools/bcftools

if [ $# -ne 2 ]; then
  echo "usage: $0 <reference indexed with bwa index> <novoindex reference>"
  exit 65
fi

if [ ! -e $1 ]; then
  echo "can't find reference .fasta: $1, please supply a bwa-indexed .fasta"
  exit 65
fi

if [ ! -e $2 ]; then
  echo "can't find novoalign index: $2, please supply a novoindex index"
  exit 65
fi

if [ ! -e $1.bwt ]; then
  echo "can't find $1.bwt: is $1 indexed with bwa?"
  exit 65
fi

python3 -O ../bin/addsv.py -p 1 -v ../test_data/test_sv.txt -f ../test_data/testregion_novo_realign.bam -r $1 -o ../test_data/testregion_sv_mut.bam -c ../test_data/test_cnvlist.txt.gz --aligner novoalign --alignopts novoref:$2
if [ $? -ne 0 ]; then
  echo "python3 -O ../bin/addsv.py failed."
  exit 65
else
  echo "sorting output bam..."
  samtools sort -T ../test_data/testregion_sv_mut.sorted.bam -o ../test_data/testregion_sv_mut.sorted.bam ../test_data/testregion_sv_mut.bam
  mv ../test_data/testregion_sv_mut.sorted.bam ../test_data/testregion_sv_mut.bam

  echo "indexing output bam..."
  samtools index ../test_data/testregion_sv_mut.bam

  echo "making pileups..."
  samtools mpileup $2 ../test_data/testregion_sv_mut.bam ../test_data/testregion_realign.bam >test_sv.pileup.txt
  echo "done. output in test_sv.pileup.txt"
fi
