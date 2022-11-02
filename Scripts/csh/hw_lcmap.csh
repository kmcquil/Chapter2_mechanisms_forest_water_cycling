#!/bin/tcsh
#BSUB -W 2:00
#BSUB -n 32
#BSUB -R span[hosts=1]
#BSUB -R "rusage[mem=20GB]"
##BSUB -x
#BSUB -o out.%J
#BSUB -e err.%J
#BSUB -J lcmap

conda activate /usr/local/usrapps/klmarti3/python-gis
python /share/klmarti3/kmcquil/Chapter2_mechanisms_forest_water_cycling/Scripts/headwater_lcmap.py
conda deactivate

