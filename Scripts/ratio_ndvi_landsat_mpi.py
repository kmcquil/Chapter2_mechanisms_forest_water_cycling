# Calculate NDVI metrics for a headwater catchment
# compatible with pynote launcher to take advantage of mpi since we have 117,000 headwater catchments 
# arg1 must be specified in the command line execution of the python script 
# arg1 is the headwater nhdplus_id 

import os
import geopandas as gpd
from glob import glob
import numpy as np
import pandas as pd 
import rasterio 
from rasterio.merge import merge
from rasterio.plot import show
from rasterio.mask import mask
import shutil
import sys

#home = "/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling"
home = "/share/klmarti3/kmcquil/Chapter2_mechanisms_forest_water_cycling"
exec(open(os.path.join(home, 'Scripts', "Functions.py")).read())

# get the hw id from the command line argument 
arg1= int(float(sys.argv[1]))

# read in the headwater catchments and subset to the catchment of interest
#hws = gpd.read_file(glob(os.path.join(home, "Data", "Catchments", "Headwater", "headwater_catchments_sbr.gpkg"))[0])
hws = gpd.read_file(os.path.join(home, "Data", "Catchments", "Headwater", "headwater_catchments_perm_forest_32617.shp"))
hw = hws[hws.NHDPlusID == arg1]

# read in the upslope contributing area raster
uca = glob(os.path.join(home, "Data", "Topography", "UCA", "uca_10m_sbr.tif"))[0]

# list all ndvi files - there are 1283
ndvi_files = glob(os.path.join(home, "Data", "NDVI","Landsat", "Landsat_images","*.tif"))

# save a csv with the metrics for each clear sky image 
catchment_ndvi_ratio_ts(hw, arg1, uca, ndvi_files)

