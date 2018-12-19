#!/bin/sh

# adds up to 100 SNPs to a ~770 kb region around the LARGE gene
# requires samtools/bcftools

REF=../test_data/Homo_sapiens_chr22_assembly19.fasta

command -v addsv.py >/dev/null 2>&1 || { echo "addsv.py isn't installed" >&2; exit 65; }

if [ ! -e $REF ]
then
    echo "can't find reference .fasta: $REF, please supply a bwa-indexed .fasta"
    exit 65
fi

if [ ! -e $REF.bwt ]
then
    echo "can't find $REF.bwt: is $REF indexed with bwa?"
    exit 65
fi

addsv.py -p 1 -v ../test_data/test_sv.txt -f ../test_data/testregion_realign.bam -r $REF -o ../test_data/testregion_sv_mut.bam --aligner mem --keepsecondary --seed 1234 --inslib ../test_data/test_inslib.fa --require_exact

if [ $? -ne 0 ]
then
  echo "addsv.py failed."
  exit 65
else
  echo "sorting output bam..."
  samtools sort -T ../test_data/testregion_sv_mut.sorted.bam -o ../test_data/testregion_sv_mut.sorted.bam ../test_data/testregion_sv_mut.bam
  mv ../test_data/testregion_sv_mut.sorted.bam ../test_data/testregion_sv_mut.bam

  echo "indexing output bam..."
  samtools index ../test_data/testregion_sv_mut.bam
fi
