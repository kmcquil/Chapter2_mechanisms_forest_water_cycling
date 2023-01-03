import os
import geopandas as gpd
from glob import glob
import numpy as np
import pandas as pd 
import rasterio 
from rasterio.merge import merge
from rasterio.plot import show
from rasterio.mask import mask
import matplotlib.pyplot as plt
import shutil
import sys

home = "/share/klmarti3/kmcquil/Chapter2_mechanisms_forest_water_cycling"
exec(open(os.path.join(home, "Scripts","Functions.py")).read())

# get the hw id from the command line argument 
arg1= int(float(sys.argv[1]))

# read in the headwater catchments and subset to the catchment of interest
hws = gpd.read_file(os.path.join(home, "Data", "Catchments", "Headwater", "headwater_catchments_perm_forest_32617.shp"))
hw = hws[hws.NHDPlusID == arg1]

# read in the upslope contributing area raster
uca = glob(os.path.join(home, "Data", "Topography", "UCA", "uca_10m_sbr.tif"))[0]

# list all of the et files - there's 150ish of them 
et_files = glob(os.path.join(home, "Data", "Evapotranspiration", "Ecostress", "Data", "clean", "*.tif"))

catchment_et_metrics(hw, arg1, uca, et_files)
