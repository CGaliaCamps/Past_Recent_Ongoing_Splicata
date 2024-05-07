#!/bin/bash

cami=/piec2/cgalia/Styela_mundial/2bRAD_world/11adaptacio/
gff=splicata_round3_lnc_sno_trna_notracks.gff
GOterm=geneid2goterm.txt

#Start the loop for multiple files
cd /piec2/cgalia/Styela_mundial/2bRAD_world/11adaptacio/

for file in *.out;

	do
		#set the root of your file names
		base=${file%%.*}

		#transfor from txt to bam file
		dos2unix $file
		sed -i 's/_/\t/2' $file
		awk -F'\t' 'BEGIN {OFS="\t"} {print $0, $2+1}' $file > $base\.bam

		#select from your gff those SNPs under selection
		/home/soft/bedtools2-2.30.0/bedtools intersect -a $cami/$gff -b $cami/$base\.bam > $cami/$base\_gene_hit.txt

		#Select the name of your genes
		cut -f 9 $cami/$base\_gene_hit.txt | cut -d ";" -f 1 | cut -d "=" -f 2 | cut -d ":" -f 1 | sort | uniq > $cami/$base\_geneID_hit.txt

		#transform from geneID to GOterm
		grep -wf $cami/$base\_geneID_hit.txt $cami/$GOterm | cut -f 2 > $base\_GOTerm_list.txt

		##in case of enrichment generate of the rest
		grep -vwf $cami/$base\_geneID_hit.txt $cami/$GOterm | cut -f 1 > $base\_geneID_NOhit.txt
	done
