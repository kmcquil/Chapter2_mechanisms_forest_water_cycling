## Use the DEM tiles from the USGS NED 1/3 arc second dataset and create a mosaiced elevation raster.
## Use the elevation to calculate aspect, slope, and TWI 

import os
import glob as glob
import numpy as np
from osgeo import gdal
import richdem as rd
import geopandas as gpd

#####################################################################################################################
# Elevation
# Create a virtual raster from the elevation tiles and then save the virtual raster to a .tif 
#home = "/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling"
home = "/share/klmarti3/kmcquil/Chapter2_mechanisms_forest_water_cycling"
elevation_tiles = glob.glob(os.path.join(home, "Data", "Topography", "Elevation", "USGS_NED_10m", "USGS*.tif" ))

vrt_options = gdal.BuildVRTOptions(resampleAlg = "bilinear", addAlpha=False)
dst_name = os.path.join(home, "Data", "Topography", "Elevation", "elevation_10m.vrt")
dem_vrt = gdal.BuildVRT(dst_name, elevation_tiles, options=vrt_options)
dem_vrt = None

dst_name_tif = os.path.join(home, "Data", "Topography", "Elevation", "elevation_10m.tif")
dem = gdal.Translate(dst_name_tif, dst_name)
dem = None

# crop the raster to the blue ridge ecoregion to reduce size
dst_name_tif_sbr = os.path.join(home, "Data", "Topography", "Elevation", "elevation_10m_sbr.tif")
sbr = os.path.join(home, "Data", "ROI", "blue_ridge_plus_reference_xl_4269.shp")
gdal.Warp(dst_name_tif_sbr, dst_name_tif,  cutlineDSName=sbr, cropToCutline=True)

dst_32617 = os.path.join(home, "Data", "Topography", "Elevation", "elevation_10m_sbr_32617.tif")
gdal.Warp(dst_32617, dst_name_tif_sbr ,dstSRS='EPSG:32617', resampleAlg = "bilinear")

#####################################################################################################################
# Aspect
dem = rd.LoadGDAL(dst_32617)
aspect = rd.TerrainAttribute(dem, attrib='aspect')
outfileaspect = os.path.join(home, "Data", "Topography", "Aspect", "aspect_10m_sbr.tif")
rd.SaveGDAL(outfileaspect, aspect)

#####################################################################################################################
# Slope
slope = rd.TerrainAttribute(dem, attrib='slope_degrees')
outfileslope = os.path.join(home, "Data", "Topography", "Slope", "slope_10m_sbr.tif")
rd.SaveGDAL(outfileslope, slope)


#####################################################################################################################
# TWI 

def twi(file): 
    # TWI =  Ln((Accumulated flow * 1 * 1)/Tan(slope radians))
    # Hwang 2012 used D-infinity method which is what we do here 
    # make sure the none of the accumulated cells are <= 0 because the Ln will mess up - realistically it should not ever be less than 1
    # and make the slope is never eactly 0 because you can't divide by 0 
    
    # first reproject to projected crs (32617)
    out_file = os.path.join(home, "Data", "Topography", "TWI", "DEM_tiles_32617", os.path.basename(file))
    gdal.Warp(out_file,file,dstSRS='EPSG:32617', resampleAlg = "bilinear")

    dem = rd.LoadGDAL(out_file)
    # fill the dem
    rd.FillDepressions(dem, epsilon=True, in_place=True)

    slope_radians = rd.TerrainAttribute(dem, attrib="slope_radians")
    slope_radians[slope_radians == slope_radians.no_data] = np.nan
    slope_radians[slope_radians <= 0] = 0.1

    flow_acc = rd.FlowAccumulation(dem, method="Dinf")
    flow_acc[flow_acc == flow_acc.no_data] = np.nan 
    flow_acc[flow_acc <= 0] = 1

    # save the flow accumulation separately for upslope contributing area 
    ucaout = os.path.join(home, "Data", "Topography", "UCA","UCA_tiles", "uca_" + os.path.basename(file))
    rd.SaveGDAL(ucaout, flow_acc)
    
    twi = np.log(flow_acc)/np.tan(slope_radians)  
    outfiletwi = os.path.join(home, "Data", "Topography", "TWI", "TWI_tiles", "twi_" + os.path.basename(file))
    rd.SaveGDAL(outfiletwi, twi)  

    # calculate TWI 
for tile in elevation_tiles:
    twi(tile)

# reproject shapefile to match the new projection of the TWI (32617)
sbr = gpd.read_file(os.path.join(home, "Data", "ROI", "blue_ridge_xl_simplified", "SBR_XL_simplified.shp"))
sbr_32617 = sbr.to_crs("EPSG:32617")
sbr_32617.to_file(os.path.join(home, "Data", "ROI", "blue_ridge_xl_simplified", "SBR_XL_simplified_32617.shp"))

# create a vrt of the twi files 
twi_tiles = glob.glob(os.path.join(home, "Data", "Topography", "TWI", "TWI_tiles", "*.tif"))
vrt_options = gdal.BuildVRTOptions(resampleAlg = "bilinear", addAlpha=False, srcNodata='nan')
dst_name = os.path.join(home, "Data", "Topography", "TWI", "twi_10m.vrt")
twi_vrt = gdal.BuildVRT(dst_name, twi_tiles, options=vrt_options)
twi_vrt = None

# crop to the cutline (sbr) and save as .tif 
dst_name_tif_sbr = os.path.join(home, "Data", "Topography", "TWI", "twi_10m_sbr.tif")
sbr = os.path.join(home, "Data", "ROI", "blue_ridge_xl_simplified", "SBR_XL_simplified_32617.shp")
gdal.Warp(dst_name_tif_sbr, dst_name,  cutlineDSName=sbr, cropToCutline=True)

# create a vrt of the uca files 
uca_tiles = glob.glob(os.path.join(home, "Data", "Topography", "UCA", "UCA_tiles", "*.tif"))
vrt_options = gdal.BuildVRTOptions(resampleAlg = "bilinear", addAlpha=False, srcNodata='nan')
dst_name = os.path.join(home, "Data", "Topography", "UCA", "uca_10m.vrt")
uca_vrt = gdal.BuildVRT(dst_name, uca_tiles, options=vrt_options)
uca_vrt = None

# crop to the cutline (sbr) and save as .tif 
dst_name_tif_sbr = os.path.join(home, "Data", "Topography", "UCA", "uca_10m_sbr.tif")
sbr = os.path.join(home, "Data", "ROI", "blue_ridge_xl_simplified", "SBR_XL_simplified_32617.shp")
gdal.Warp(dst_name_tif_sbr, dst_name,  cutlineDSName=sbr, cropToCutline=True)