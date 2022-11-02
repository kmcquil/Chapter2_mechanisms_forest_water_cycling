bsub -Is -n 1 -W 30 tcsh
hostname   #This should not be a login node!
module load conda
conda activate /usr/local/usrapps/klmarti3/python-gis
python