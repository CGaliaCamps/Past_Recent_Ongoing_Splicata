#qsub -q "*@node100" -l "h_vmem=100G" p_value.sh  ###to submit the work
#!/bin/bash

export PATH=/home/soft/gdal-2.3.2/bin/:$PATH

export LD_LIBRARY_PATH=/home/soft/gdal-2.3.2/lib/:$LD_LIBRARY_PATH

/home/soft/R-4.1.0/bin/R --vanilla </home/usuaris/cgalia/.Scripts/Styela_2bRAD/permuting_FST_pvalue_NCSC.R > /piec2/cgalia/Styela_mundial/2bRAD_world/pvalue_FST/FST_pvalue.out 2 > /piec2/cgalia/Styela_mundial/2bRAD_world/pvalue_FST/FST_pvalue.err


