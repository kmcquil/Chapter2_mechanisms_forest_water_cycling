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


shp_file = os.path.join(home, "Data", "Catchments", "Headwater", "headwater_catchments_sbr_shp_4269.shp")
shp = gpd.read_file(os.path.join(home, "Data", "Catchments", "Headwater", "headwater_catchments_sbr_shp_4269.shp"))
#shp_file = os.path.join(home, "Data", "Catchments", "Headwater", "TESTING_headwater_catchments_sbr_shp_4269.shp")
#shp = gpd.read_file(os.path.join(home, "Data", "Catchments", "Headwater", "TESTING_headwater_catchments_sbr_shp_4269.shp"))
with fiona.open(shp_file) as src:
        features = list(src)
print("shapes done!")

tifs = [os.path.join(home, "Data", "Topography", "Elevation", "elevation_10m_sbr.tif"),
        os.path.join(home, "Data", "Topography", "Slope", "slope_10m_sbr.tif"),
        os.path.join(home, "Data", "Topography", "Aspect", "aspect_10m_sbr.tif"), 
        os.path.join(home, "Data", "latitude.tif")]
colnames = ["Elevation", "Slope", "Aspect", "Latitude"]

K = -1
if __name__ == "__main__":
    
    # Create a process pool using all cores
    cores = 36
    p = multiprocessing.Pool(cores)

    for tif in tifs:
        print(tif)
        tic = time.perf_counter()
        K +=1
        stats_lists = p.map(partial(zonal_stats_partial, tif), chunks(features, cores))

        # flatten to a single list
        stats = pd.DataFrame(list(itertools.chain(*stats_lists)))
        stats.columns = [colnames[K]]
        assert len(stats) == len(features)

        df = pd.concat([shp['NHDPlusID'], stats], axis=1)
        df.to_csv(os.path.join(home, "Data", "Topography", colnames[K] + "_hw_summary.csv"))
        toc = time.perf_counter()
        print('iteration ' + str(K) + ' ran for ' + str(toc-tic) + ' seconds')
    
    p.close()
    p.join()
