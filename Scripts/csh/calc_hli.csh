#!/bin/tcsh
#BSUB -W 12:00
#BSUB -n 24
#BSUB -R span[hosts=1]
#BSUB -R "rusage[mem=50GB]"
##BSUB -x 
#BSUB -o out.%J
#BSUB -e err.%J
#BSUB -J hli

module load R/4.0.2-gcc4.8.5
setenv LD_LIBRARY_PATH /usr/local/apps/gdal/gcc483-2.1.1/lib:$LD_LIBRARY_PATH 
setenv PATH /usr/local/apps/gdal/gcc483-2.1.1/bin:$PATH 
setenv LD_LIBRARY_PATH /usr/local/apps/proj/gcc483-4.9.3/lib:$LD_LIBRARY_PATH 
Rscript /share/klmarti3/kmcquil/Chapter2_mechanisms_forest_water_cycling/Scripts/calculate_hli.R
