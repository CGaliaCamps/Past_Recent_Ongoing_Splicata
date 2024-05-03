#!/bin/bash

#########################     FOR READS PROPERLY TRIMMED    ######################

GENOME="/piec2/cgalia/GenomaSplicata/Splicata_final_assembly_smasked.fasta"
READS="/piec2/cgalia/Styela_mundial/2bRAD_world/2fastaxbRAD/Trimmed"
OUTPUT="/piec2/cgalia/Styela_mundial/2bRAD_world/3genotyping/Definitius/smasked"


###index your genome###
bwa index ${GENOME}


cd $READS

mkdir ${OUTPUT%/*}
mkdir $OUTPUT


for file in *.fasta;

	do
		name=${file%.*}

		###map, sort and save as bam your alignment###

		bwa mem ${GENOME} ${READS}/$file -t 18 -M  > ${OUTPUT}/$name\_map.sam
		/home/soft/samtools-1.3.1/samtools addreplacerg ${OUTPUT}/$name\_map.sam -r ID:$name -r SM:$name > ${OUTPUT}/$name\_rg.sam
		rm ${OUTPUT}/$name\_map.sam
		/home/soft/samtools-1.3.1/samtools view -S -b ${OUTPUT}/$name\_rg.sam > ${OUTPUT}/$name\.bam
		rm ${OUTPUT}/$name\_rg.sam
		/home/soft/samtools-1.3.1/samtools sort ${OUTPUT}/$name\.bam > ${OUTPUT}/$name\_sorted.bam
		/home/soft/samtools-1.3.1/samtools index ${OUTPUT}/$name\_sorted.bam

	done



