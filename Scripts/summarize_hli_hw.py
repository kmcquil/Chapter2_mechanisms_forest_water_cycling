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
    return zonal_stats(feats, tif, all_touched=True, stats='mean')

shp_file = os.path.join(home, "Data", "Catchments", "Headwater", "headwater_catchments_perm_forest_32617.shp")
shp = gpd.read_file(shp_file)
with fiona.open(shp_file) as src:
        features = list(src)
print("shapes done!")

tif = os.path.join(home, "Data", "Topography", "HLI_max.tif")

if __name__ == "__main__":
    
    # Create a process pool using all cores
    cores = 24
    p = multiprocessing.Pool(cores)
    stats_lists = p.map(partial(zonal_stats_partial, tif), chunks(features, cores))
    # flatten to a single list
    stats = pd.DataFrame(list(itertools.chain(*stats_lists)))
    df = pd.concat([shp['NHDPlusID'], stats], axis=1)
    df.to_csv(os.path.join(home, "Data", "Topography", "HLI_hw_summary.csv"))

