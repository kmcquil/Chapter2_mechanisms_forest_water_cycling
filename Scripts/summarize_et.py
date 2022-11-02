# Summarize average ET metrics from all of the clear sky scenes 

import os
from glob import glob
import numpy as np
import pandas as pd 
import shutil
import sys

home = "/share/klmarti3/kmcquil/Chapter2_mechanisms_forest_water_cycling"
#home = "/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling"
exec(open(os.path.join(home, 'Scripts', "Functions.py")).read())
files = glob(os.path.join(home, "Data", "Evapotranspiration", "catchment_et_metrics", "*.csv"))

row_list = []
for file in files:
    print(file)
    results = summarize_et_catchment(file)
    if results == None:
        continue
    row_list.append(results)
summary_df = pd.DataFrame(row_list)
summary_df.to_csv(os.path.join(home, "Data","Evapotranspiration","summary_et.csv"))



