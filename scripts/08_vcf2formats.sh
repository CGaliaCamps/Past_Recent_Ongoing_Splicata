#!/bin/bash

cd /piec2/cgalia/Styela_mundial/2bRAD_world/5vcf2formats

mkdir structure
mkdir plink
mkdir nexus
mkdir fasta
mkdir phylip
mkdir binary_nexus
mkdir treemix

source vcf2phylip/bin/activate
pip install numpy

for file in *.vcf;

	do

		base="${file%%.recode.vcf}"

####get your treemix popmap

		bcftools query -l $file > IDs_$base.txt
                sed 's/\(..\).*/\1/' IDs_$base.txt > pops_$base.txt
		paste IDs_$base.txt pops_$base.txt | column -s $'\t' -t > treemix_popmap_$base.txt


#####get your converted files

		#to structure
		/piec1/software/plink-1.9/plink --allow-extra-chr --const-fid --vcf $file --recode structure --out ./structure/$base.str
		#to phylip
		python /piec1/software/vcf2phylip/vcf2phylip.py -i $file -r --output-folder ./phylip
		#to nexus
		python /piec1/software/vcf2phylip/vcf2phylip.py -i $file -r -n -p --output-folder ./nexus
		#to fasta
		python /piec1/software/vcf2phylip/vcf2phylip.py -i $file -r -p -f --output-folder ./fasta
		#to plink
		/home/soft/vcftools/bin/vcftools --vcf $file --plink --out ./plink/$base.plink
		#to binary nexus
		python /piec1/software/vcf2phylip/vcf2phylip.py	-i $file -r -p -b --output-folder ./binary_nexus
		#to treemix
                python /piec1/software/vcf2treemix.py $file treemix_popmap_$base.txt ./treemix/$base.treemix
		#to 012

	done

cd plink

for file in *.ped;

	do
		base=${file%.ped}
		nom=${file%.recode.plink}

		/piec1/software/plink-1.9/plink --noweb --file $base --recodeA --out $nom

	done


