# Create headwater dataset using NHDPlus HR 

import fiona 
import geopandas as gpd
import numpy as np
import pandas as pd
import glob as glob
from osgeo import gdal

# List the .gdb files 
home = "/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling"
gdb_files = glob.glob(os.path.join(home, "Data","Catchments", "Headwater", "NHDPlus_HR", "*", "*.gdb"))

# Get all the layers from the .gdb file to figure out which ones we need to identify headwaters  
layers = fiona.listlayers(gdb_files[0])

# function to identify and save headwaters within the roi 
def identify_headwaters(file):
    # read in the catchment and properties associated with flow lines 
    catchment = gpd.read_file(file,layer="NHDPlusCatchment")
    flowlineVAA = gpd.read_file(file, layer = "NHDPlusFlowlineVAA")
    flowlineVAA = flowlineVAA.drop('geometry',1)

    # subset to just headwater flowlines 
    hw_flowlineVAA = flowlineVAA[flowlineVAA['HWType'] == 0.0]

    # merge with the catchments and only keep catchments that have a hw flowline 
    hw_catchments = catchment.merge(hw_flowlineVAA, on="NHDPlusID", how = 'inner') 

    # keep the catchment within the roi (roi being the blue ridge ecoregion + the refreence watersheds that extend slightly )
    roi = gpd.read_file(os.path.join(home, "Data", "ROI", "blue_ridge_plus_reference.shp"))
    roi = roi.to_crs(hw_catchments.crs)
    hw_catchments_join= gpd.tools.sjoin(hw_catchments, roi, predicate = "intersects", how = "left")
    hw_catchments_keep = hw_catchments_join.dropna(subset=['index_right'])

    cols = ['NHDPlusID', 'AreaSqKm_x', 'StreamLeve', 'StreamOrde','geometry']
    keep = hw_catchments_keep[cols] 

    if keep.shape[0] == 0: 
        return 0
    else: 
        # save the result
        keep.to_file(os.path.join(home, "Data","Catchments", "Headwater", "Headwaters_ROI", os.path.basename(file)[:-4] + ".gpkg"), driver="GPKG")
        return keep.shape[0]

# loop through each gdb file to create headwaters for each section (HUC4)
total_size = []
for file in gdb_files:
    total_size.append(identify_headwaters(file))
    print('Finished: ' + file)

# create one gpkg with all headwaters 
hw_files = glob(os.path.join(home, "Data", "Catchments", "Headwater", "Headwaters_ROI", "*.gpkg"))
i = 0
for file in hw_files: 
    i = i+1
    if i == 1: 
        hws = gpd.read_file(file)
    else:
        hw = gpd.read_file(file)
        hws = pd.concat([hws, hw])

# repoject to match topography data 
hws = hws.to_crs('EPSG:32617')

# save it 
hws.to_file(os.path.join(home, "Data", "Catchments", "Headwater", "headwater_catchments_sbr.gpkg"))