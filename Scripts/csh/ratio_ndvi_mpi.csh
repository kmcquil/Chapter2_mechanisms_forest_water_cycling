#!/bin/tcsh
#BSUB -n 100  # number of MPI processes
#BSUB -W 72:00  # maximum time
#BSUB -R "rusage[mem=1GB]"
#BSUB -oo tasks_out
#BSUB -eo tasks_err
#BSUB -J tasks  # job name

module load PrgEnv-intel
module load conda
conda activate /usr/local/usrapps/klmarti3/python-gis
mpiexec python -m mpi4py -m pynodelauncher /share/klmarti3/kmcquil/Chapter2_mechanisms_forest_water_cycling/Scripts/commands_ratio_ndvi_landsat_mpi_permanent_forest.txt
