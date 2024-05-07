#!/bin/bash

PATH=/piec2/cgalia/Styela_mundial/2bRAD_world/3genotyping/max-minDP/

cd $PATH

mkdir filtered
mkdir filtered/regions

for file in *.vcf;
	do
		base=${file%%.vcf}

#filter by abundance
		/home/soft/vcftools/bin/vcftools --vcf $PATH/$file --remove-indels --minGQ 98 --minDP 5 --min-alleles 2 --max-alleles 2 --mac 2 --recode --out $PATH/filtered/2mac_5DP_$base

#get the mean depth statistics
		/home/soft/vcftools/bin/vcftools --vcf $PATH/filtered/2mac_5DP_$base.recode.vcf --site-mean-depth --out $PATH/filtered/2mac_5DP_$base
	done

#filter by depth (please check which heading referes to depth in your case)
/home/soft/vcftools/bin/vcftools --vcf $PATH/filtered/2mac_5DP_2bRAD_nou_raw.recode.vcf --max-meanDP 172 --recode --out $PATH/filtered/2mac_5-172DP_2bRAD_nou

#filter by presence
/home/soft/vcftools/bin/vcftools --vcf $PATH/filtered/2mac_5-172DP_2bRAD_nou.recode.vcf --max-missing 0.7 --recode --out $PATH/filtered/2mac_5-172DP_70_2bRAD_nou

cd filtered/

for file in *5-172DP*;
	do
		base=${file%%.recode.vcf}

		/home/soft/vcftools/bin/vcftools --vcf $PATH/filtered/$file --remove-indv NC14 --remove-indv NC18 --remove-indv NC16 --remove-indv NC15 --remove-indv NC25 --recode --out ./regions/$base\_NoNC
		/home/soft/vcftools/bin/vcftools --vcf $PATH/filtered/$file --remove-indv NC14 --remove-indv NC18 --remove-indv NC16 --remove-indv NC15 --remove-indv NC25 --remove-indv SC4 --remove-indv SC11 --remove-indv SC24 --remove-indv SC18 --remove-indv SC23 --recode --out ./regions/$base\_NoNCSC
		/home/soft/vcftools/bin/vcftools --vcf $PATH/filtered/$file --keep atl_names.txt --recode --out ./regions/$base\_atl
		/home/soft/vcftools/bin/vcftools --vcf $PATH/filtered/$file --keep pac_names.txt --recode --out ./regions/$base\_pac
	done


cd /piec2/cgalia/Styela_mundial/2bRAD_world/3genotyping/max-minDP/filtered/regions/

base=Chromosome_

mkdir Chr_No_inv

for file in *vcf;
	do
		nom=${file%%.recode.vcf}
		/home/soft/vcftools/bin/vcftools --vcf $file --chr $base\1 --chr $base\3 --chr $base\5 --chr $base\6 --chr $base\7 --chr $base\8 --chr $base\9 --chr $base\10 --chr $base\12 --chr $base\13 --chr $base\14 --chr $base\15 --recode --out ./Chr_No_inv/noInv_$nom
	done


for file in *.vcf;
	do
		for n in 2 4 11 16;
			do
				mkdir chr$n
				nom=${file%%.recode.vcf}
				/home/soft/vcftools/bin/vcftools --vcf $file --chr $base\2 --recode --out ./chr2/2\_$nom
				/home/soft/vcftools/bin/vcftools --vcf $file --chr $base\4 --recode --out ./chr4/4\_$nom
				/home/soft/vcftools/bin/vcftools --vcf $file --chr $base\11 --recode --out ./chr11/11\_$nom
				/home/soft/vcftools/bin/vcftools --vcf $file --chr $base\16 --recode --out ./chr16/16\_$nom
			done
	done

