import os
import time
from glob import glob

import rasterio as rio
from osgeo import gdal, osr
from rasterstats import zonal_stats
import fiona
import geopandas as gpd
import pandas as pd
import numpy as np

import multiprocessing
from functools import partial
import itertools

home = "/share/klmarti3/kmcquil/Chapter2_mechanisms_forest_water_cycling"
#home = "/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling"

def chunks(data, n):
    """Yield successive n-sized chunks from a slice-able iterable."""
    for i in range(0, len(data), n):
        yield data[i:i+n]

def zonal_stats_partial(tif, feats):
    """Wrapper for zonal stats, takes a list of features"""
    return zonal_stats(feats, tif, all_touched=True, stats='mean')

def get_file_info(files):
    row_list = []
    for file in files:
        filename = file
        basename = os.path.basename(filename)
        breakdown_path = basename.split("_")   

        Type = breakdown_path[0]
        Var = breakdown_path[1]
        if Type == "anomalies": 
            Season = breakdown_path[2]
            Year = breakdown_path[3][:-4]
        else: 
            Season = breakdown_path[2][:-4]
            Year = np.nan 

        row_list.append({'Type':Type, 'Var':Var, 'Season':Season, 'Year':Year})
    
    return pd.DataFrame(row_list)

def summarize_multithread(pool, variable):
    
    files = glob(os.path.join(home, "Data", "Climate", "tifs", "*" + variable  + "*.tif"))
    file_info = get_file_info(files)
    print(len(files))
    K = -1

    for tif in files:
        
        # track progress 
        print(tif)
        tic = time.perf_counter()
        K +=1
        
        # map the function to chunks of features across the requested cores and convert to a df 
        stats_lists = pool.map(partial(zonal_stats_partial, tif), chunks(features, cores))
        stats = pd.DataFrame(list(itertools.chain(*stats_lists)))
        stats.columns = ['Value']
        df = pd.concat([shp['NHDPlusID'], stats], axis=1)

        # add additional info for indexing 
        df['Type'] = [file_info.iloc[K].Type]* df.shape[0]
        df['Var'] = [file_info.iloc[K].Var]* df.shape[0]
        df['Season'] = [file_info.iloc[K].Season]* df.shape[0]
        df['Year'] = [file_info.iloc[K].Year]* df.shape[0]

        if K == 0:
            df_all = df
        else:
            df_all = pd.concat([df_all, df], axis=0)
        
        toc = time.perf_counter()
        print('iteration ' + str(K) + ' ran for ' + str(toc-tic) + ' seconds')

    p.close()
    p.join()
    df_all.to_csv(os.path.join(home, "Data", "Climate", "Summary", variable + ".csv"))


shp_file = os.path.join(home, "Data", "Catchments", "Headwater", "headwater_catchments_perm_forest_climate.shp")
shp = gpd.read_file(shp_file)

with fiona.open(shp_file) as src:
        features = list(src)
print("shapes done!")

if __name__ == "__main__":
    # Create a process pool using all cores
    cores = 36
    p = multiprocessing.Pool(cores)
    summarize_multithread(p, 'drydays')