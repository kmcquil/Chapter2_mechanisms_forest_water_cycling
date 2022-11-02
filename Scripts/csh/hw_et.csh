#!/bin/tcsh
#BSUB -W 02:00
#BSUB -n 1
#BSUB -R span[hosts=1]
#BSUB -R "rusage[mem=10GB]"
##BSUB -x 
#BSUB -o out.%J
#BSUB -e err.%J
#BSUB -J et

conda activate /usr/local/usrapps/klmarti3/python-gis
python /share/klmarti3/kmcquil/Chapter2_mechanisms_forest_water_cycling/Scripts/summarize_et.py
conda deactivate