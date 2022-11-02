from rasterstats import zonal_stats
import fiona
import geopandas as gpd
import os
import pandas as pd
import itertools
import multiprocessing

home = "/share/klmarti3/kmcquil/Chapter2_mechanisms_forest_water_cycling"

# reproject the geopackage hw catchment file to be a shp with a crs of 4326
#gpkg = gpd.read_file(os.path.join(home, "Data", "Catchments", "Headwater", "headwater_catchments_sbr.gpkg"))
#shp = gpkg.to_crs("EPSG:4326")
#shp.to_file(os.path.join(home, "Data", "Catchments", "Headwater", "headwater_catchments_sbr_shp.shp"))

# GLOBALS
shp = os.path.join(home, "Data", "Catchments", "Headwater", "headwater_catchments_sbr_shp.shp")
tif = os.path.join(home, "Data", "lcmap_permanent_forest", "lcmap_permanent_forest.tif")


# https://github.com/perrygeo/python-rasterstats/blob/master/examples/multiproc.py
def chunks(data, n):
    """Yield successive n-sized chunks from a slice-able iterable."""
    for i in range(0, len(data), n):
        yield data[i:i+n]


def zonal_stats_partial(feats):
    """Wrapper for zonal stats, takes a list of features"""
   #return zonal_stats(feats, tif, all_touched=True, stats = ['min'])
    return zonal_stats(feats, tif, all_touched=True, categorical=True)

if __name__ == "__main__":
    
    with fiona.open(shp) as src:
        features = list(src)

    # Create a process pool using all cores
    cores = 32
    p = multiprocessing.Pool(cores)

    # parallel map
    stats_lists = p.map(zonal_stats_partial, chunks(features, cores))

    # flatten to a single list
    stats = list(itertools.chain(*stats_lists))

    assert len(stats) == len(features)


# put back into a dataframe, calculate the pct forest, and save to a csv
df = pd.DataFrame(stats)
#df['percent_forest'] = df.iloc[:,1] / df.sum(axis=1)
shape_gdf = gpd.read_file(shp)["NHDPlusID"]
df_concat = pd.concat([df, shape_gdf], axis=1)
df_concat.to_csv(os.path.join(home, "Data", "lcmap_permanent_forest", "lcmap_headwater_pct_forest.csv"))



# also find the permanent forest composition for all of the reference watersheds 
shp = os.path.join(home, "Data", "Catchments", "Reference", "gages_ii",  "reference_keep_4326.shp")
with fiona.open(shp) as src:
    features = list(src)
tif = os.path.join(home, "Data", "lcmap_permanent_forest", "lcmap_permanent_forest.tif")

stats_lists = zonal_stats(features, tif, all_touched=True, categorical=True)
df = pd.DataFrame(stats_lists)
shape_gdf = gpd.read_file(shp)["GAGE_ID"]
df_concat = pd.concat([df, shape_gdf], axis=1)
df_concat.to_csv(os.path.join(home, "Data", "lcmap_permanent_forest", "lcmap_reference_watersheds.csv"))

