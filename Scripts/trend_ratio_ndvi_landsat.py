import os
from glob import glob
import numpy as np
import pandas as pd 
import shutil
import sys
#import multiprocess
#import itertools

home = "/share/klmarti3/kmcquil/Chapter2_mechanisms_forest_water_cycling"
#home = "/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling"
exec(open(os.path.join(home, 'Scripts', "Functions.py")).read())
files = glob(os.path.join(home, "Data", "NDVI", "catchment_ratio_ndvi_results", "*.csv"))

row_list = []
for file in files:
    print(file)
    results = trend_ratio_ndvi(file)
    if results == None:
        continue
    row_list.append(results[0])
    row_list.append(results[1])
    row_list.append(results[2])
    row_list.append(results[3])
trends_df = pd.DataFrame(row_list)
trends_df.to_csv(os.path.join(home, "Data", "NDVI", "ratio_ndvi_trend_results.csv"))

