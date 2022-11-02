#!/bin/tcsh
#BSUB -n 100  # number of MPI processes
#BSUB -W 72:00  # maximum time
#BSUB -R "rusage[mem=1GB]"
#BSUB -oo tasks_out_et
#BSUB -eo tasks_err_et
#BSUB -J tasks_et  # job name

module load PrgEnv-intel
module load conda
conda activate /usr/local/usrapps/klmarti3/python-gis
mpiexec python -m mpi4py -m pynodelauncher /share/klmarti3/kmcquil/Chapter2_mechanisms_forest_water_cycling/Scripts/commands_et_metrics_mpi_permanent_forest.txt 

