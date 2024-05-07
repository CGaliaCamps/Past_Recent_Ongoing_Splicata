#!/bin/bash


PATH=/piec2/cgalia/Styela_mundial/2bRAD_world/max-meanDP

cd $PATH

for folder in */plink/;
        do
                cd $folder
                mkdir MDS

                for file in *.map;
                        do
                                file2=${file%%.map}
				base=${file2%%.plink}

				/piec1/software/plink-1.9/plink --noweb --file $file2  --make-bed --out preMDS_$base
				/piec1/software/plink-1.9/plink --noweb --bfile preMDS_$base --cluster --mds-plot 50 --out ./MDS/MDS_$base
                        done

		cd $PATH
        done


PATH=/piec2/cgalia/Styela_mundial/2bRAD_world/12neutralstructure/vcf2plink/plink

cd $PATH

                for file in *.map;
                        do
                                file2=${file%%.map}
				base=${file2%%.plink}

				/piec1/software/plink-1.9/plink --noweb --file $file2  --make-bed --out preMDS_$base
				/piec1/software/plink-1.9/plink --noweb --bfile preMDS_$base --cluster --mds-plot 50 --out ./MDS/MDS_$base
                        done





