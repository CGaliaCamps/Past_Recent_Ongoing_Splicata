#!/bin/bash

export PATH=/usr/lib/jvm/java-9-openjdk-amd64/bin/:$PATH
export LD_LIBRARY_PATH=/usr/lib/jvm/java-9-openjdk-amd64/lib/:$LD_LIBRARY_PATH

INPATH=/piec1/cgalia/Styela_mundial/WGS_medit/1rawdata/
OUTPATH=/piec1/cgalia/Styela_mundial/WGS_medit/1filtdata/

gunzip ${INPATH}/*.gz

cd ${INPATH}

for file in  *_1.fq;
do
	name=${file%_*}
	java -jar /home/soft/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 16 -phred33 ${INPATH}/$name\_1.fq ${INPATH}/$name\_2.fq ${OUTPATH}/$name\_forward_paired.fq ${OUTPATH}/$name\_forward_unpaired.fq.gz ${OUTPATH}/$name\_reverse_paired.fq ${OUTPATH}/$name\_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

gzip ${INPATH}/*.fq


