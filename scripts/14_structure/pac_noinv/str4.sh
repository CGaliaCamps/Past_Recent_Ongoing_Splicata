#! /bin/bash
path=/piec2/cgalia/Styela_mundial/2bRAD_world/7structure/pac_noinv_nomissing/

cd $path

for file in *str;

	do
		base=${file%%.str}
		number=$(sed 1d $file | cut -f 1 | sort | uniq | wc -l);
		loci=$(head -n 1 $file | awk '{FS=" "} ; {print NF}');

		for t in {1..10};

			do

#			/piec1/software/Structure/console/structure -m /piec1/software/Structure/console/mainparams -e /piec1/software/Structure/console/extraparams -i $path$file -D 123 -L $loci -N $number -K 1 -o $path/$file\_1_$t.txt > $path/$file\_K1_$t &
#			/piec1/software/Structure/console/structure -m /piec1/software/Structure/console/mainparams -e /piec1/software/Structure/console/extraparams -i $path/$file -D 123 -L $loci -N $number -K 2 -o $$path/$base\_2_$t.txt > $path/$base\_K2_$t &
#			/piec1/software/Structure/console/structure -m /piec1/software/Structure/console/mainparams -e /piec1/software/Structure/console/extraparams -i $path$file -D 123 -L $loci -N $number -K 3 -o $path/$file\_3_$t.txt > $path/$file\_K3_$t &
			/piec1/software/Structure/console/structure -m /piec1/software/Structure/console/mainparams -e /piec1/software/Structure/console/extraparams -i $path$file -D 123 -L $loci -N $number -K 4 -o $path/$file\_4_$t.txt > $path/$file\_K4_$t &
#			/piec1/software/Structure/console/structure -m /piec1/software/Structure/console/mainparams -e /piec1/software/Structure/console/extraparams -i $path$file -D 123 -L $loci -N $number -K 5 -o $path/$file\_5_$t.txt > $path/$file\_K5_$t &
#			/piec1/software/Structure/console/structure -m /piec1/software/Structure/console/mainparams -e /piec1/software/Structure/console/extraparams -i $path$file -D 123 -L $loci -N $number -K 6 -o $path/$file\_6_$t.txt > $path/$file\_K6_$t &
#			/piec1/software/Structure/console/structure -m /piec1/software/Structure/console/mainparams -e /piec1/software/Structure/console/extraparams -i $path$file -D 123 -L $loci -N $number -K 7 -o $path/$file\_7_$t.txt > $path/$file\_K7_$t &
#			/piec1/software/Structure/console/structure -m /piec1/software/Structure/console/mainparams -e /piec1/software/Structure/console/extraparams -i $path$file -D 123 -L $loci -N $number -K 8 -o $path/$file\_8_$t.txt > $path/$file\_K8_$t &
#			/piec1/software/Structure/console/structure -m /piec1/software/Structure/console/mainparams -e /piec1/software/Structure/console/extraparams -i $path$file -D 123 -L $loci -N $number -K 9 -o $path/$file\_9_$t.txt > $path/$file\_K9_$t &
#			/piec1/software/Structure/console/structure -m /piec1/software/Structure/console/mainparams -e /piec1/software/Structure/console/extraparams -i $path$file -D 123 -L $loci -N $number -K 10 -o $path/$file\_10_$t.txt > $path/$file\_K10_$t & 
#			/piec1/software/Structure/console/structure -m /piec1/software/Structure/console/mainparams -e /piec1/software/Structure/console/extraparams -i $path$file -D 123 -L $loci -N $number -K 11 -o $path/$file\_11_$t.txt > $path/$file\_K11_$t 
#			/piec1/software/Structure/console/structure -m /piec1/software/Structure/console/mainparams -e /piec1/software/Structure/console/extraparams -i $path$file -D 123 -L $loci -N $number -K 12 -o $path/$file\_12_$t.txt > $path/$file\_K12_$t 
#			/piec1/software/Structure/console/structure -m /piec1/software/Structure/console/mainparams -e /piec1/software/Structure/console/extraparams -i $path$file -D 123 -L $loci -N $number -K 13 -o $path/$file\_13_$t.txt > $path/$file\_K13_$t 
#			/piec1/software/Structure/console/structure -m /piec1/software/Structure/console/mainparams -e /piec1/software/Structure/console/extraparams -i $path$file -D 123 -L $loci -N $number -K 14 -o $path/$file\_14_$t.txt > $path/$file\_K14_$t 
#			/piec1/software/Structure/console/structure -m /piec1/software/Structure/console/mainparams -e /piec1/software/Structure/console/extraparams -i $path$file -D 123 -L $loci -N $number -K 15 -o $path/$file\_15_$t.txt > $path/$file\_K15_$t 
#			/piec1/software/Structure/console/structure -m /piec1/software/Structure/console/mainparams -e /piec1/software/Structure/console/extraparams -i $path$file -D 123 -L $loci -N $number -K 16 -o $path/$file\_16_$t.txt > $path/$file\_K16_$t 
#			/piec1/software/Structure/console/structure -m /piec1/software/Structure/console/mainparams -e /piec1/software/Structure/console/extraparams -i $path$file -D 123 -L $loci -N $number -K 17 -o $path/$file\_17_$t.txt > $path/$file\_K17_$t 
#			/piec1/software/Structure/console/structure -m /piec1/software/Structure/console/mainparams -e /piec1/software/Structure/console/extraparams -i $path$file -D 123 -L $loci -N $number -K 18 -o $path/$file\_18_$t.txt > $path/$file\_K18_$t 
#			/piec1/software/Structure/console/structure -m /piec1/software/Structure/console/mainparams -e /piec1/software/Structure/console/extraparams -i $path$file -D 123 -L $loci -N $number -K 19 -o $path/$file\_19_$t.txt > $path/$file\_K19_$t 
#			/piec1/software/Structure/console/structure -m /piec1/software/Structure/console/mainparams -e /piec1/software/Structure/console/extraparams -i $path$file -D 123 -L $loci -N $number -K 20 -o $path/$file\_20_$t.txt > $path/$file\_K20_$t 

			wait

			done
	done

