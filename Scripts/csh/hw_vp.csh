#!/bin/tcsh
#BSUB -W 4:00
#BSUB -n 24
#BSUB -R span[hosts=1]
#BSUB -R "rusage[mem=30GB]"
##BSUB -x
#BSUB -o out.%J
#BSUB -e err.%J
#BSUB -J vp

conda activate /usr/local/usrapps/klmarti3/python-gis
python /share/klmarti3/kmcquil/Chapter2_mechanisms_forest_water_cycling/Scripts/summarize_gs_vp_hw.py
conda deactivate
