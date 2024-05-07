#!/bin/bash

cd /piec2/cgalia/Styela_mundial/2bRAD_world/8demography/get_neutral/sfs_stairway

python /piec1/software/easySFS.py -i neutral_loci_all.vcf -p global_popmap.sfs --preview -a 1> global_preview.txt
python /piec1/software/easySFS.py -i neutral_loci_all.vcf -p region_popmap.sfs --preview -a 1> region_preview.txt
python /piec1/software/easySFS.py -i neutral_loci_all.vcf -p single_popmap.sfs --preview -a 1> single_preview.txt

python /piec1/software/easySFS.py -i neutral_loci_all.vcf -p global_popmap.sfs -a --proj 122 -f -o global_sfs_max
python /piec1/software/easySFS.py -i neutral_loci_all.vcf -p region_popmap.sfs -a --proj 48,48,12 -f -o regions_sfs_max
python /piec1/software/easySFS.py -i neutral_loci_all.vcf -p single_popmap.sfs -a --proj 2,8,6,8,8,8,4,8,8,10,6,8,8,8,8,8,8,8 -f -o single_sfs_max

python /piec1/software/easySFS.py -i neutral_loci_all.vcf -p global_popmap.sfs -a --proj 152 -f -o global_sfs_trade
python /piec1/software/easySFS.py -i neutral_loci_all.vcf -p region_popmap.sfs -a --proj 63,60,16 -f -o regions_sfs_trade
python /piec1/software/easySFS.py -i neutral_loci_all.vcf -p single_popmap.sfs -a --proj 2,10,8,10,10,10,6,10,10,12,8,10,10,10,8,8,10,10 -f -o single_sfs_trade
