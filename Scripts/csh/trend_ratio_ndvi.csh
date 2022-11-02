#!/bin/tcsh
#BSUB -W 02:00
#BSUB -n 16
#BSUB -R span[hosts=1]
#BSUB -R "rusage[mem=30GB]"
##BSUB -x 
#BSUB -o out.%J
#BSUB -e err.%J

conda activate /usr/local/usrapps/klmarti3/python-gis
python /share/klmarti3/kmcquil/Chapter2_mechanisms_forest_water_cycling/Scripts/trend_ratio_ndvi_landsat.py
conda deactivate