#! /bin/bash
path=/piec2/cgalia/Styela_mundial/2bRAD_world/5vcf2formats/structure/

cd $path

for file in *str;

	do
		base=${file%%.str}
		number=$(sed 1d $file | cut -f 1 | sort | uniq | wc -l);
		loci=$(head -n 1 $file | awk '{FS=" "} ; {print NF}');

		for t in {1..10};

			do

#			/piec1/software/Structure/console/structure -m /piec1/software/Structure/console/mainparams -e /piec1/software/Structure/console/extraparams -i $path$file -D 123 -L $loci -N $number -K 1 -o /piec2/cgalia/Styela_mundial/2bRAD_world/7structure/$file\_1_$t.txt > /piec2/cgalia/Styela_mundial/2bRAD_world/7structure/$file\_K1_$t &
#			/piec1/software/Structure/console/structure -m /piec1/software/Structure/console/mainparams -e /piec1/software/Structure/console/extraparams -i $path/$file -D 123 -L $loci -N $number -K 2 -o $/piec2/cgalia/Styela_mundial/2bRAD_world/7structure/$base\_2_$t.txt > /piec2/cgalia/Styela_mundial/2bRAD_world/7structure/$base\_K2_$t &
#			/piec1/software/Structure/console/structure -m /piec1/software/Structure/console/mainparams -e /piec1/software/Structure/console/extraparams -i $path$file -D 123 -L $loci -N $number -K 3 -o /piec2/cgalia/Styela_mundial/2bRAD_world/7structure/$file\_3_$t.txt > /piec2/cgalia/Styela_mundial/2bRAD_world/7structure/$file\_K3_$t &
#			/piec1/software/Structure/console/structure -m /piec1/software/Structure/console/mainparams -e /piec1/software/Structure/console/extraparams -i $path$file -D 123 -L $loci -N $number -K 4 -o /piec2/cgalia/Styela_mundial/2bRAD_world/7structure/$file\_4_$t.txt > /piec2/cgalia/Styela_mundial/2bRAD_world/7structure/$file\_K4_$t &
#			/piec1/software/Structure/console/structure -m /piec1/software/Structure/console/mainparams -e /piec1/software/Structure/console/extraparams -i $path$file -D 123 -L $loci -N $number -K 5 -o /piec2/cgalia/Styela_mundial/2bRAD_world/7structure/$file\_5_$t.txt > /piec2/cgalia/Styela_mundial/2bRAD_world/7structure/$file\_K5_$t &
#			/piec1/software/Structure/console/structure -m /piec1/software/Structure/console/mainparams -e /piec1/software/Structure/console/extraparams -i $path$file -D 123 -L $loci -N $number -K 6 -o /piec2/cgalia/Styela_mundial/2bRAD_world/7structure/$file\_6_$t.txt > /piec2/cgalia/Styela_mundial/2bRAD_world/7structure/$file\_K6_$t &
#			/piec1/software/Structure/console/structure -m /piec1/software/Structure/console/mainparams -e /piec1/software/Structure/console/extraparams -i $path$file -D 123 -L $loci -N $number -K 7 -o /piec2/cgalia/Styela_mundial/2bRAD_world/7structure/$file\_7_$t.txt > /piec2/cgalia/Styela_mundial/2bRAD_world/7structure/$file\_K7_$t &
#			/piec1/software/Structure/console/structure -m /piec1/software/Structure/console/mainparams -e /piec1/software/Structure/console/extraparams -i $path$file -D 123 -L $loci -N $number -K 8 -o /piec2/cgalia/Styela_mundial/2bRAD_world/7structure/$file\_8_$t.txt > /piec2/cgalia/Styela_mundial/2bRAD_world/7structure/$file\_K8_$t &
#			/piec1/software/Structure/console/structure -m /piec1/software/Structure/console/mainparams -e /piec1/software/Structure/console/extraparams -i $path$file -D 123 -L $loci -N $number -K 9 -o /piec2/cgalia/Styela_mundial/2bRAD_world/7structure/$file\_9_$t.txt > /piec2/cgalia/Styela_mundial/2bRAD_world/7structure/$file\_K9_$t &
#			/piec1/software/Structure/console/structure -m /piec1/software/Structure/console/mainparams -e /piec1/software/Structure/console/extraparams -i $path$file -D 123 -L $loci -N $number -K 10 -o /piec2/cgalia/Styela_mundial/2bRAD_world/7structure/$file\_10_$t.txt > /piec2/cgalia/Styela_mundial/2bRAD_world/7structure/$file\_K10_$t & 
#			/piec1/software/Structure/console/structure -m /piec1/software/Structure/console/mainparams -e /piec1/software/Structure/console/extraparams -i $path$file -D 123 -L $loci -N $number -K 11 -o /piec2/cgalia/Styela_mundial/2bRAD_world/7structure/$file\_11_$t.txt > /piec2/cgalia/Styela_mundial/2bRAD_world/7structure/$file\_K11_$t 
#			/piec1/software/Structure/console/structure -m /piec1/software/Structure/console/mainparams -e /piec1/software/Structure/console/extraparams -i $path$file -D 123 -L $loci -N $number -K 12 -o /piec2/cgalia/Styela_mundial/2bRAD_world/7structure/$file\_12_$t.txt > /piec2/cgalia/Styela_mundial/2bRAD_world/7structure/$file\_K12_$t 
#			/piec1/software/Structure/console/structure -m /piec1/software/Structure/console/mainparams -e /piec1/software/Structure/console/extraparams -i $path$file -D 123 -L $loci -N $number -K 13 -o /piec2/cgalia/Styela_mundial/2bRAD_world/7structure/$file\_13_$t.txt > /piec2/cgalia/Styela_mundial/2bRAD_world/7structure/$file\_K13_$t 
#			/piec1/software/Structure/console/structure -m /piec1/software/Structure/console/mainparams -e /piec1/software/Structure/console/extraparams -i $path$file -D 123 -L $loci -N $number -K 14 -o /piec2/cgalia/Styela_mundial/2bRAD_world/7structure/$file\_14_$t.txt > /piec2/cgalia/Styela_mundial/2bRAD_world/7structure/$file\_K14_$t 
#			/piec1/software/Structure/console/structure -m /piec1/software/Structure/console/mainparams -e /piec1/software/Structure/console/extraparams -i $path$file -D 123 -L $loci -N $number -K 15 -o /piec2/cgalia/Styela_mundial/2bRAD_world/7structure/$file\_15_$t.txt > /piec2/cgalia/Styela_mundial/2bRAD_world/7structure/$file\_K15_$t 
#			/piec1/software/Structure/console/structure -m /piec1/software/Structure/console/mainparams -e /piec1/software/Structure/console/extraparams -i $path$file -D 123 -L $loci -N $number -K 16 -o /piec2/cgalia/Styela_mundial/2bRAD_world/7structure/$file\_16_$t.txt > /piec2/cgalia/Styela_mundial/2bRAD_world/7structure/$file\_K16_$t 
#			/piec1/software/Structure/console/structure -m /piec1/software/Structure/console/mainparams -e /piec1/software/Structure/console/extraparams -i $path$file -D 123 -L $loci -N $number -K 17 -o /piec2/cgalia/Styela_mundial/2bRAD_world/7structure/$file\_17_$t.txt > /piec2/cgalia/Styela_mundial/2bRAD_world/7structure/$file\_K17_$t 
#			/piec1/software/Structure/console/structure -m /piec1/software/Structure/console/mainparams -e /piec1/software/Structure/console/extraparams -i $path$file -D 123 -L $loci -N $number -K 18 -o /piec2/cgalia/Styela_mundial/2bRAD_world/7structure/$file\_18_$t.txt > /piec2/cgalia/Styela_mundial/2bRAD_world/7structure/$file\_K18_$t 
			/piec1/software/Structure/console/structure -m /piec1/software/Structure/console/mainparams -e /piec1/software/Structure/console/extraparams -i $path$file -D 123 -L $loci -N $number -K 19 -o /piec2/cgalia/Styela_mundial/2bRAD_world/7structure/$file\_19_$t.txt > /piec2/cgalia/Styela_mundial/2bRAD_world/7structure/$file\_K19_$t 
#			/piec1/software/Structure/console/structure -m /piec1/software/Structure/console/mainparams -e /piec1/software/Structure/console/extraparams -i $path$file -D 123 -L $loci -N $number -K 20 -o /piec2/cgalia/Styela_mundial/2bRAD_world/7structure/$file\_20_$t.txt > /piec2/cgalia/Styela_mundial/2bRAD_world/7structure/$file\_K20_$t 

			wait

			done
	done

