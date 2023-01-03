import os
from glob import glob
import numpy as np
import pandas as pd 
import rasterio 
from rasterio.merge import merge
from rasterio.plot import show
from rasterio.mask import mask
import shutil
import sys
import multiprocess
import itertools

home = "/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling"
exec(open(os.path.join(home, "Scripts","Functions.py")).read())

# list the qc flags and et files to clean  
qcdf = pd.read_csv(os.path.join(home, "Data", "Evapotranspiration", "Ecostress", "Meta", "ECO2LSTE-001-SDS-QC-lookup.csv"))
flags = qcdf[qcdf["Mandatory QA flags"] == "Pixel produced, best quality"].Value
et_files = glob(os.path.join(home, "Data", "Evapotranspiration", "Ecostress", "Data", "*ETinst*.tif"))

# define the function to clean each file and save 
def clean_et(file):

    scene_id = file[len(file)-28:len(file)-12]  
    qc_scene = glob(os.path.join(home, "Data", "Evapotranspiration", "Ecostress", "Data","ECO2LSTE"+"*"+scene_id+"*"+".tif"))[0] 
    qc = rasterio.open(qc_scene).read(1)
    pixel_flag_mask = np.isin(qc, flags)

    with rasterio.open(file, GDAL_DISABLE_READDIR_ON_OPEN=False) as src:
        img = np.ma.masked_values(src.read(masked=True), np.nan)
        img = np.where(pixel_flag_mask, img, np.nan)

        out_name = os.path.join(home, "Data", "Evapotranspiration", "Ecostress", "Data", "clean", os.path.basename(file))
        out_meta = src.meta
        out_meta.update({"driver": "GTiff",
                 "nodata": -9999})
        with rasterio.open(out_name, "w", **out_meta) as dest:
            dest.write(img)


if __name__ == "__main__":
       
    # Create a process pool using all cores
    cores = multiprocess.cpu_count() - 2
    p = multiprocess.Pool(cores)
    print(p)

    # parallel map
    p.map(clean_et, et_files)


