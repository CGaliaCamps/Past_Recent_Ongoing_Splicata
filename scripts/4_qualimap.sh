#!/bin/bash

export LD_LIBRARY_PATH=/home/soft/R-4.1.0/lib/:$LD_LIBRARY_PATH
export PATH=/home/soft/R-4.1.0/bin/:$PATH

path=/piec1/cgalia/Styela_mundial/2bRAD_world/3genotyping/bwamapped/

cd $path

for file in *sorted.bam;

	do
		base="${file%%_*}"

		/home/soft/qualimap/qualimap_v2.2.1/qualimap bamqc -bam $path/$file -nt 4 -outdir $path/qualimaps/$base/ -outfile $base\_map -outformat HTML -c
	done

