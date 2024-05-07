#!/bin/bash

source /home/soft/virtenvs/python3/multiqc/.mqcvenv/bin/activate
export LANG=C.UTF-8


#TEST raw data

multiqc /piec2/cgalia/Styela_mundial/2bRAD_world/1raw_data/fastQC/* -o /piec2/cgalia/Styela_mundial/2bRAD_world/1raw_data/fastQC/multiqc/


#TEST mappings#

INPATH=/piec1/cgalia/Styela_mundial/2bRAD_world/3genotyping/bwamapped/qualimaps/
OUTPATH=/piec1/cgalia/Styela_mundial/2bRAD_world/3genotyping/bwamapped/qualimaps/TOTSqualimaps/

multiqc $INPATH/*/* -o $OUTPATH
