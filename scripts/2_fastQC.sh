#!/bin/bash

RAW=/piec1/cgalia/Styela_mundial/WGS_medit/1rawdata/
FILT=/piec1/cgalia/Styela_mundial/WGS_medit/1filtdata/

mkdir $RAW/fastQC/
mkdir $FILT/fastQC/

/home/soft/FastQC/fastqc $RAW/*.fq -t 12 -o $RAW/fastQC/
/home/soft/FastQC/fastqc $FILT/15/*.fq -t 12 -o $FILT/fastQC/

