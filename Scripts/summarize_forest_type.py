from rasterstats import zonal_stats
import rasterio
import fiona
import geopandas as gpd
import os
import pandas as pd
import itertools
import multiprocessing
from functools import partial
import time

home = "/share/klmarti3/kmcquil/Chapter2_mechanisms_forest_water_cycling"

def chunks(data, n):
    """Yield successive n-sized chunks from a slice-able iterable."""
    for i in range(0, len(data), n):
        yield data[i:i+n]

def zonal_stats_partial(tif, feats):
    """Wrapper for zonal stats, takes a list of features"""
    return zonal_stats(feats, tif, all_touched=True, categorical=True)

shp_file = os.path.join(home, "Data", "Catchments", "Headwater", "headwater_catchments_perm_forest.shp")
shp = gpd.read_file(shp_file)
with fiona.open(shp_file) as src:
        features = list(src)
print("shapes done!")

tifs = [os.path.join(home, "Data", "nlcd_permanent_forest", "nlcd_2001.tif"),
        os.path.join(home, "Data", "nlcd_permanent_forest", "nlcd_2016.tif")]
colnames = ["nlcd_2001", "nlcd_2016"]

if __name__ == "__main__":
    
    # Create a process pool using all cores
    cores = 24
    p = multiprocessing.Pool(cores)

    tif = tifs[0]
    nam = colnames[0]
    stats_lists = p.map(partial(zonal_stats_partial, tif), chunks(features, cores))
    # flatten to a single list
    stats = pd.DataFrame(list(itertools.chain(*stats_lists)))
    df = pd.concat([shp['NHDPlusID'], stats], axis=1)
    df.to_csv(os.path.join(home, "Data", "nlcd_permanent_forest", "hw_" + nam +"_summary.csv"))

    tif = tifs[1]
    nam = colnames[1]
    stats_lists = p.map(partial(zonal_stats_partial, tif), chunks(features, cores))
    # flatten to a single list
    stats = pd.DataFrame(list(itertools.chain(*stats_lists)))
    df = pd.concat([shp['NHDPlusID'], stats], axis=1)
    df.to_csv(os.path.join(home, "Data", "nlcd_permanent_forest", "hw_" + nam +"_summary.csv"))
    
    p.close()
    p.join()
