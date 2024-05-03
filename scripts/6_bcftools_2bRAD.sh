#!/bin/bash

INPATH=/piec2/cgalia/Styela_mundial/2bRAD_world/3genotyping/max-minDP/
OUTPATH=/piec2/cgalia/Styela_mundial/2bRAD_world/3genotyping/max-minDP/
REF=/piec2/cgalia/GenomaSplicata/Splicata_final_assembly.fasta
CORES=8

mkdir $OUTPATH

realpath $INPATH/*sorted.bam  > $OUTPATH/pathlist.txt

/home/soft/samtools-20211022/bcftools/bcftools mpileup --threads $CORES -Ov --annotate FORMAT/AD,FORMAT/DP,INFO/AD --fasta-ref $REF -b $OUTPATH/pathlist.txt | /home/soft/samtools-20211022/bcftools/bcftools call --threads $CORES -mv | /home/soft/samtools-20211022/bcftools/bcftools annotate --set-id +'%CHROM\_%POS' > ${OUTPATH}/2bRAD_raw.vcf
