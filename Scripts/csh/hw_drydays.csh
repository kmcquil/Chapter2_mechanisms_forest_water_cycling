#!/bin/tcsh
#BSUB -W 24:00
#BSUB -n 36
#BSUB -R span[hosts=1]
#BSUB -R "rusage[mem=30GB]"
##BSUB -x
#BSUB -o out.%J
#BSUB -e err.%J
#BSUB -J drydays

conda activate /usr/local/usrapps/klmarti3/python-gis
python /share/klmarti3/kmcquil/Chapter2_mechanisms_forest_water_cycling/Scripts/summarize_drydays_hw.py
conda deactivate
