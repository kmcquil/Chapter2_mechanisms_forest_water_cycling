# Function names are all lower case and have _ between words 
# Variable names are also all lower case and have _ between words 
# Constants should be capitalized with _ between words 
# Classes should be CamelCase 

#  Collection of functions used in the ET + HVG analysis 
import os
import re
import geopandas as gpd
from glob import glob
import numpy as np
import pandas as pd 
import rasterio 
from rasterio.merge import merge
from rasterio.plot import show
from rasterio.mask import mask
import richdem as rd 
from osgeo import gdal
import fiona
from sklearn.linear_model import LinearRegression
from scipy import stats
import shutil
import time

def check_missing_data(image, boundary):
    # Check if there is data missing in an image within a shapefile boundary 
    # Inputs 
    #   image = filepath of raster 
    #   boundary = geometry of boundary - usually from geopanadas object but idk if it matters 
    # Outputs
    #   If there is no missing data, return the cropped image and metadata
    #   Otherwise, return False

    # open image and boundary and make sure they completely overlap 
    image = rasterio.open(image, GDAL_DISABLE_READDIR_ON_OPEN=True)
    img_extent = str(image.bounds)
    #p = re.compile(r'\d+\.\d+')  # Compile a pattern to capture float values
    p = re.compile(r'[+-]?\d+\.\d+') 
    tif_minx, tif_miny, tif_maxx, tif_maxy = [float(i) for i in p.findall(img_extent)]  # Convert strings to float
    # make sure the boundary is in the same crs as the image 
    boundary=boundary.to_crs(image.crs)
    shp_minx, shp_miny, shp_maxx, shp_maxy = boundary.geometry.total_bounds

    if ( (shp_minx < tif_minx) | (shp_miny < tif_miny) | (shp_maxx > tif_maxx) | (shp_maxy > tif_maxy)):
        return False
    else: 
        # crop the image to the boundary 
        out_image, out_transform = rasterio.mask.mask(image, boundary, crop=True, nodata= -9999)
        # find how many empy pixels are within the boundary 
        empties = np.sum(np.isnan(out_image))
        if(empties == 0):
            return out_image, out_transform
        else:
            return False 

def get_epsg(img):
    # Retrieve the EPSG code from a raster
    # Input
    #   img: a raster opened with rasterio (eg img = rasterio.open(image))
    # Output 
    #   code: a string with the EPSG code 
    crs = img.crs
    if crs.is_epsg_code:
        code = int(crs['init'].lstrip('epsg:'))
    return code

def check_metadata(image, original, out_image, out_transform):
    # Check if the crs, extent, and resolution of two images match 
    # Input
    #   imgage: a raster file path 
    #   template: a raster file path for comparison
    # Output 
    #   logical of True or False  

    img = rasterio.open(image, GDAL_DISABLE_READDIR_ON_OPEN=True)
    og = rasterio.open(original, GDAL_DISABLE_READDIR_ON_OPEN=True)
    # check CRS
    img_crs = get_epsg(img)
    temp_crs = get_epsg(og)
    crs_match = img_crs == temp_crs

    # check resolution 
    img_res = img.res
    temp_res = og.res
    res_match = img_res == temp_res

    # check extent 
    img_extent = str(img.bounds)
    p = re.compile(r'\d+\.\d+')  # Compile a pattern to capture float values
    floats = [float(i) for i in p.findall(img_extent)]  # Convert strings to float
    
    xres = out_transform[0]
    yres = -out_transform[4]
    ulx = out_transform[2]
    uly = out_transform[5]
    xsize = np.shape(out_image)[1]
    ysize = np.shape(out_image)[2]
    lrx = ulx + (ysize * xres)
    lry = uly - (xsize * yres)
    temp_extent = [ulx, lry, lrx, uly]
    extent_match = floats == temp_extent

    if ((crs_match == True) and (extent_match == True) and (res_match == True)):
        return True
    else: 
        return False


def resample_to_match(image, outfile, template, alg):
    # Resample an image to match a tmemplate 
    # Inputs: 
    #   image = image file path to be resampled 
    #   outfile = where to write the resmapled image 
    #   template = image filepath that serves as a template that we want the scene to match (extent, resolution, crs)
    #   alg = algorithm to use for resampling
    # Ouputs: 
    #   write out the new resampled image to the outfile 
    reference = gdal.Open(template)  # this opens the file in only reading mode
    band = reference.GetRasterBand(1)
    
    ulx, xres, xskew, uly, yskew, yres  = reference.GetGeoTransform() # get the extent and resolution 
    yres=-yres # make sure yres is positive 
    lrx = ulx + (reference.RasterXSize * xres)
    lry = uly - (reference.RasterYSize * yres)
    outputbounds =  [ulx, lry, lrx, uly] # output bounds as (minX, minY, maxX, maxY) in target SRS
    dest = band.GetNoDataValue()
    dstSRS= reference.GetProjection()
    src = gdal.Open(image)
    s_band = src.GetRasterBand(1)
    sest = s_band.GetNoDataValue()
    
    kwargs = {"format": "GTiff", "xRes":xres, "yRes":yres, "outputBounds":outputbounds, 
        "resampleAlg":alg, "dstNodata":dest, "srcNodata":sest, "dstSRS":dstSRS}
    gdal.Warp(outfile, image, **kwargs)


def crop_to_boundary(image, boundary, outname):
    # Crop and mask an image to a polygon
    # Inputs: 
    #   scene = filepath of the scene you want to crop 
    #   boundary = shapefile of the boundary you want to crop to 
    #   outname = where to save the new cropped image 
    # Outputs: 
    #   Write out the cropped raster 

    # with fiona.open(boundary, "r") as shapefile:
    # shapes = [feature["geometry"] for feature in shapefile]
    shapes = boundary['geometry']
    with rasterio.open(image, GDAL_DISABLE_READDIR_ON_OPEN=True) as src:
        out_image, out_transform = rasterio.mask.mask(src, shapes, crop=True)
        out_meta = src.meta
        outName=outname
        out_meta.update({"driver": "GTiff",
                 "height": out_image.shape[1],
                 "width": out_image.shape[2],
                 "transform": out_transform})
        with rasterio.open(outName, "w", **out_meta, GDAL_DISABLE_READDIR_ON_OPEN=True) as dest:
            dest.write(out_image)

# similar to crop_to_boundary but takes the filename as the boundary input  
def crop_to_shp(image, boundary, outname):
    # Crop and mask an image to a polygon 
    # Inputs: 
    #   scene = filepath of the scene to crop 
    #   boundary = filepath of the shapefile to crop to the boundary 
    #   outname = filepath to save the new cropped image 
    with fiona.open(boundary, "r") as shapefile:
        shapes = [feature["geometry"] for feature in shapefile]
    with rasterio.open(image, GDAL_DISABLE_READDIR_ON_OPEN=True) as src:
        out_image, out_transform = rasterio.mask.mask(src, shapes, crop=True)
        out_meta = src.meta
        outName=outname
        out_meta.update({"driver": "GTiff",
                 "height": out_image.shape[1],
                 "width": out_image.shape[2],
                 "transform": out_transform})
        with rasterio.open(outName, "w", **out_meta, GDAL_DISABLE_READDIR_ON_OPEN=True) as dest:
            dest.write(out_image)

def summarize_ndvi_twi(ndvi, twi):    # formerly vhg
    # Create a dataframe of the average, sd, and count of NDVI in binned TWI 
    # Inputs: 
    #   ndvi = np.array of ndvi 
    #   twi = np.array of twi matching ndvi 
    # Outputs: 
    #   results = dataframe of the average, sd, count of NDVI, TWI bin
    ndvi = rasterio.open(ndvi, GDAL_DISABLE_READDIR_ON_OPEN=True).read(1)
    mask = ndvi == 0
    ndvi = np.ma.array(ndvi, mask=mask)

    twi = rasterio.open(twi, GDAL_DISABLE_READDIR_ON_OPEN=True).read(1)
    twi = np.ma.array(twi, mask=mask)

    # chop the twi tail at the 95 quantile 
    q95 = np.quantile(twi.compressed(), 0.9)
    def round_off(number): 
        return round(number*2)/2
    q95 = round_off(q95)
    length = int((2*(q95-2)) + 1)
    start_interval = list(np.linspace(2, q95, length))

    ndvi_avg_interval = []
    ndvi_sd_interval = []
    ndvi_count_interval = []
    
    for i in start_interval: 
        interval_mask = ((twi < i) | (twi >= (i+0.5))) # this is what we want to mask out --- so the opposite of the interval we are interested in 
        ndvi_interval = np.ma.array(ndvi, mask=interval_mask)
        ndvi_count_interval.append(ndvi_interval.count(axis=None))
        if ndvi_interval.count(axis=None) < 15:
            ndvi_avg_interval.append(np.nan)
            ndvi_sd_interval.append(np.nan)
        else:
            ndvi_avg_interval.append(np.nanmean(ndvi_interval, axis=None))
            ndvi_sd_interval.append(np.nanstd(ndvi_interval, axis=None))
        
    end_interval = list(map(lambda x : x + 0.5, start_interval))
    results = {
        'Start_TWI':start_interval, 
        'End_TWI':end_interval,
        'Binned_NDVI_mean':ndvi_avg_interval,
        'Binned_NDVI_sd':ndvi_sd_interval,
        'Binned_NDVI_count':ndvi_count_interval
    }
    results = pd.DataFrame(results)
    return results 


def calculate_catchment_ndvi_twi_summary(ndvi, ndvi_out, twi, twi_out, headwater, headwaterID):   # formerly vhgheadwater
    # Steps: 
    #   Crop and mask NDVI to a headwater catchment
    #   Resample, crop, and mask TWI to match the NDVI for the headwater catchment 
    #   Calculate metrics of the binned NDVI and return as a dataframe
    # Inputs 
    #   ndvi = filepath of ndvi 
    #   ndvi_out = filepath to write the ndvi that is cropped to a headwater 
    #   twi = filepath of twi 
    #   twi_out = filepath to write the twi that is resampled/cropped to the headwater ndvi 
    #   headwater = boundary of the catchment to target 
    #   headwaterID = string of headwater id 
    # Outputs
    #   
    crop_to_boundary(ndvi, headwater, ndvi_out)
    resample_to_match(twi, twi_out, ndvi_out, "average")
    vhg = summarize_ndvi_twi(ndvi_out, twi_out)
    vhg["WSID"] = [headwaterID] * vhg.shape[0]
    return vhg



def calc_topography(dempath, boundary, landscape, WSID):
    # check the dem and shapefile are both in the same projection. reproject if necessary. 
    # put a 100m buffer around the shapefile
    # crop the dem to the buffered shapefile 
    # calculate aspect and slope on the cropped dem and save 
    # calculate TWI on the cropped dem and save 

    # read in the dem and ensure the boundary matches the crs 
    # they need to be in a projected coordinate system for buffer calculation 
    dem=rasterio.open(dempath, GDAL_DISABLE_READDIR_ON_OPEN=True)
    sub=boundary.to_crs(dem.crs)

    # make a new folder for topo data for this headwater and save the boundary to the new folder
    os.makedirs(os.path.join(os.getcwd(), "Data", "Topography", landscape, WSID), exist_ok=True) # make the output directory 
    outfile_boundary = os.path.join(os.getcwd(), "Data", "Topography", landscape, WSID, "boundary.shp") # write out the file 
    sub.to_file(outfile_boundary)

    # create a buffer for later resampling/twi calcs and save 
    sub_buffer = sub.buffer(100, 16)
    outfile_buffer = os.path.join(os.getcwd(), "Data", "Topography", landscape, WSID, "boundary_buffer.shp") # write out the file 
    sub_buffer.to_file(outfile_buffer)  

    # crop the dem to the buffered boundary and save 
    outfiledem = os.path.join(os.getcwd(), "Data", "Topography", landscape, WSID, "hydem_buff.tif") # write out the file
    crop_to_shp(dempath, outfile_buffer, outfiledem)

    # calculate aspect on cropped DEM and save
    dem = rd.LoadGDAL(outfiledem)
    aspect = rd.TerrainAttribute(dem, attrib='aspect')
    outfileaspect = os.path.join(os.getcwd(), "Data", "Topography", landscape, WSID, "aspect_buff.tif") # write out the file
    rd.SaveGDAL(outfileaspect, aspect)

    # calculate TWI 
    # TWI =  Ln((Accumulated flow * 1 * 1)/Tan(slope radians))
    # Hwang 2012 used D-infinity method which is what we do here 
    # make sure the none of the accumulated cells are <= 0 because the Ln will mess up 
    slope_radians = rd.TerrainAttribute(dem, attrib="slope_radians")
    slope_radians[slope_radians == slope_radians.no_data] = np.nan
    flow_acc = rd.FlowAccumulation(dem, method="Dinf")
    flow_acc[flow_acc == flow_acc.no_data] = np.nan 
    flow_acc[flow_acc <= 0] = 0.1
    twi = np.log(flow_acc)/np.tan(slope_radians) 

    outfilesr = os.path.join(os.getcwd(), "Data", "Topography", landscape, WSID, "sloperadians_buff.tif") # write out the file
    rd.SaveGDAL(outfilesr, slope_radians)

    outfileacc = os.path.join(os.getcwd(), "Data", "Topography", landscape, WSID, "flowacc_buff.tif") # write out the file
    rd.SaveGDAL(outfileacc, flow_acc)

    outfiletwi = os.path.join(os.getcwd(), "Data", "Topography", landscape, WSID, "TWI_buff.tif") # write out the file
    rd.SaveGDAL(outfiletwi, twi)

    #rd.rdShow(twi, axes=False, cmap='magma', figsize=(8, 5.5))
    #plt.show()


def clean_ecostress(landscape, scene, flags):  # formerly cleaning
    #We used the Ecostress Level 2 LSTE Mandatory QC product. All pixels not flagged as "Pixel produced, best quality" were masked out.
    # inputs
    #   scene=the daily ET scene
    #   flags=the QC values that correspond to pixels we want to drop

    # find the mandatory qc scene that corresponds to the et scene 
    scene_id = scene[len(scene)-28:len(scene)-12]  
    qc_scene = glob.glob(os.path.join(os.getcwd(), "Data", "ETInst",landscape,"ECO2LSTE"+"*"+scene_id+"*"+".tif"))[0] 
    qc = rasterio.open(qc_scene, GDAL_DISABLE_READDIR_ON_OPEN=True) 
    qc = qc.read(1)
    qc_bool = np.isin(qc, flags)

    #mask the et scene by the qc scene. only keep pixels that match the qc flags 
    et = rasterio.open(scene, GDAL_DISABLE_READDIR_ON_OPEN=True)
    et = et.read(1)
    #et_mm = (et*0.0864)/2.45   # convert from watts/m2 to mm/day before masking 
    et_final = np.where(qc_bool, et, np.nan)


    # put back into a geotiff format and save to a new folder named clean
    # to write out - start by getting parameters
    src_dataset = gdal.Open(scene)
    src_data = src_dataset.ReadAsArray()
    geotransform = src_dataset.GetGeoTransform()
    spatialreference = src_dataset.GetProjection()
    ncol = src_dataset.RasterXSize
    nrow = src_dataset.RasterYSize
    nband = 1

    # create dataset for output
    outfile = os.path.join(os.getcwd(), "Data", "ETInst",landscape, "Clean", os.path.basename(scene))
    fmt = 'GTiff'
    driver = gdal.GetDriverByName(fmt)
    dst_dataset = driver.Create(outfile, ncol, nrow, nband,  gdal.GDT_Float32)  ##################
    dst_dataset.SetGeoTransform(geotransform)
    dst_dataset.SetProjection(spatialreference)
    # dst_dataset.GetRasterBand(1).SetNoDataValue(-9999)
    dst_dataset.GetRasterBand(1).WriteArray(et_final)  ###################
    dst_dataset = None

def clean_ecostress_catchment(landscape, boundary): # formerly cleanwhere
    # get list of daily ET files
    et_files = glob.glob(os.path.join(os.getcwd(), "Data", "ETInst",landscape,"*ETinst*.tif"))

    # apply the function to clean and convert to mm to each scene and do qc
    for x in et_files:
        clean_ecostress(landscape,x, keep)
    
    # now take the clean files and resample them to the same grid 
    clean_files = glob.glob(os.path.join(os.getcwd(), "Data", "ETInst", landscape, "Clean", "*_ETinst_*.tif"))

    # identify the parameters 
    # reference = gdal.Open(clean_files[0], 0)  # this opens the file in only reading mode
    reference = gdal.Open(clean_files[0])  # this opens the file in only reading mode
    ulx, xres, xskew, uly, yskew, yres  = reference.GetGeoTransform() # get the extent and resolution 
    yres=-yres # make sure yres is positive 
    lrx = ulx + (reference.RasterXSize * xres)
    lry = uly - (reference.RasterYSize * yres)
    outputbounds =  [ulx, lry, lrx, uly] # output bounds as (minX, minY, maxX, maxY) in target SRS
    dest = -9999999999999
    kwargs = {"format": "GTiff", "xRes":xres, "yRes":yres, "outputBounds":outputbounds, "resampleAlg":"near", "dstNodata":dest}

    # loop through each file, resample to common grid, and crop to coweeta 
    for x in clean_files:
        inputFile = x
        outputFile = os.path.join(os.getcwd(), "Data", "ETInst", landscape,"Clean", "Resampled", os.path.basename(x))
        gdal.Warp(outputFile, inputFile, **kwargs)
        crop_to_shp(outputFile, boundary, outputFile)


def calculate_uca(dem, outpath):
    # Calculate the upslope contributing area 
    # Inputs 
    #   dem = filepath of dem 
    #   outpath = folder to save the uca  
    dem = rd.LoadGDAL(dem)
    rd.FillDepressions(dem, epsilon=True, in_place=True)
    accum_dinf = rd.FlowAccumulation(dem, method='Dinf')
    out = os.path.join(outpath, "uca.tif")
    rd.SaveGDAL(out, accum_dinf)

def resample_to_match_io(image, outfile, original, out_image, out_transform, alg):
    # Resample one raster to match another raster.
    # This is a special case where the raster serving as the template has been cropped using rasterIO 
    # Inputs: 
    #   image = image to be resampled 
    #   outfile = filename for the new resampled raster 
    #   original = image serving as template BEFORE it was cropped
    #   out_image = cropped image in np.array form 
    #   out_transform = second part of the np array that comes from cropping in raterio
    # Outputs: 
    #   write out the new raster to the outfile 
    xres = out_transform[0]
    yres = -out_transform[4]

    ulx = out_transform[2]
    uly = out_transform[5]

    xsize = np.shape(out_image)[1]
    ysize = np.shape(out_image)[2]

    lrx = ulx + (ysize * xres)
    lry = uly - (xsize * yres)

    outputbounds = [ulx, lry, lrx, uly]
    
    reference = gdal.Open(original)  # this opens the file in only reading mode
    band = reference.GetRasterBand(1)
    dest = band.GetNoDataValue()
    dstSRS= reference.GetProjection()

    src = gdal.Open(image)
    s_band = src.GetRasterBand(1)
    sest = s_band.GetNoDataValue()

    kwargs = {"format": "GTiff", "xRes":xres, "yRes":yres, "outputBounds":outputbounds, 
        "resampleAlg":alg, "dstNodata":dest, "srcNodata":sest, "dstSRS":dstSRS}
    gdal.Warp(outfile, image, **kwargs)


def classify_uca(uca_filepath, ndvi_array):
    # create the upslpoe and downslope data masks 
    # Inputs: 
    # uca_filepath = uca image filepath 
    # ndvi_array = the array of ndvi values that we will use for the mask 
    # Outputs: 
    # upslope and downslope masks
    uca = rasterio.open(uca_filepath, GDAL_DISABLE_READDIR_ON_OPEN=True).read(1)
    # mask out no data values in the ndvi values to make sure we're not including unused pixels
    ndvi = ndvi_array
    mask = ndvi == -9999
    uca_match = np.ma.array(uca, mask=mask)

    # use the 75th percentile to designate upslope vs downslope 
    q75 = np.nanquantile(uca_match.compressed(), 0.75) # this has all -9999 vlaues masked so it won't count htose in the calculation
    upslope_mask = ((ndvi == -9999) | (uca >= q75))  # mask is the values you don't want 
    downslope_mask =  ((ndvi == -9999) | (uca < q75))
    return upslope_mask, downslope_mask

def calculate_ratio_ndvi(ndvi_array, upslope_mask, downslope_mask):
    # Calculate the average upslope ndvi, downslope ndvi, and ratio of downslope/upslope
    ndvi_downslope = np.ma.array(ndvi_array, mask = downslope_mask)
    ndvi_upslope = np.ma.array(ndvi_array, mask = upslope_mask)
    ndvi_downslope_avg = np.nanmean(ndvi_downslope)
    ndvi_upslope_avg = np.nanmean(ndvi_upslope)
    ratio_ndvi = ndvi_downslope_avg/ndvi_upslope_avg

    sd_mask = ndvi_array == -9999
    ndvi = np.ma.array(ndvi_array, mask = sd_mask)
    sd_ndvi = np.std(ndvi)
    mean_ndvi = np.mean(ndvi)

    results = [ndvi_upslope_avg, ndvi_downslope_avg, ratio_ndvi, sd_ndvi, mean_ndvi]
    return results

def calculate_ratio_et(et_array, upslope_mask, downslope_mask):
    # Calculate the average upslope ndvi, downslope ndvi, and ratio of downslope/upslope
    et_downslope = np.ma.array(et_array, mask = downslope_mask)
    et_upslope = np.ma.array(et_array, mask = upslope_mask)
    et_downslope_avg = np.nanmean(et_downslope)
    et_upslope_avg = np.nanmean(et_upslope)
    ratio_et = et_downslope_avg/et_upslope_avg

    sd_mask = et_array == -9999
    et = np.ma.array(et_array, mask = sd_mask)
    sd_et = np.std(et)
    avg_et = np.mean(et)

    results = [et_upslope_avg, et_downslope_avg, ratio_et, sd_et, avg_et]
    return results

def catchment_ndvi_ratio_ts(boundary, wsid, uca, ndvi_file_list):
    # Summarize upslope and downslope ndvi + the ratio for all images with no missing data and save to a csv for the given catchment
    # Inputs 
    #   boundary= catchment boundary 
    #   wsid = catchment ID 
    #   dem = DEM in the original resolution cropped + buffer around the catchment 
    #   ndvi_file_list = list of all potential ndvi files to use 
    # Outputs 
    #   a csv is saved with the average upslope, downslope, and ratio downslope:upslope ndvi for every image with no missing data for a given catchment 
    
    # create uca folder with the wsid so that we can do this for multiple watersheds at a time 
    uca_path = os.path.join(home, "Data", "Topography", "UCA","UCA_NDVI" + str(wsid))
    is_exist = os.path.exists(uca_path)
    if not is_exist:
        os.makedirs(uca_path)
    # crop the UCA to boundary and save it to the watershed specific uca folder 
    uca_path_name = os.path.join(uca_path, "uca.tif")
    crop_to_boundary(uca, boundary, uca_path_name)

    # start a list to collect rows in form of dictionaries
    row_list = []

    for image in ndvi_file_list: 
        # check for missing data within the boundary and if there is missing data skip to the next image
        # otherwise the function returns the image and transform 
        missing_data = check_missing_data(image, boundary['geometry'])
        if missing_data == False:
            continue

        # if there is no missing data, check if there is a UCA raster that matches the new cropped NDVI 
        # if not, resample the UCA to match the cropped NDVI 
        uca_files = glob(os.path.join(uca_path, "*.tif"))
        for file in uca_files:
            check = check_metadata(file, image, missing_data[0], missing_data[1])
            if check == True: 
                matching_uca = file
                break
            else:
                matching_uca = False
        if matching_uca == False:
            uca_out = os.path.join(uca_path, "uca" + str(len(uca_files) + 1) + ".tif")
            resample_to_match_io(os.path.join(uca_path, "uca.tif"), uca_out, image, missing_data[0], missing_data[1], "bilinear")
            matching_uca = uca_out
           
        # classify upslope vs downslope pixels in the resampled UCA 
        upslope_mask, downslope_mask = classify_uca(matching_uca, missing_data[0])
        
        # take the average ndvi for upslope and downslope and find the ratio 
        ndvi = calculate_ratio_ndvi(missing_data[0], upslope_mask, downslope_mask)

        # get the date 
        date = os.path.basename(image)[12:20]
        landsat = os.path.basename(image)[0:4]
        # create a dictionary 
        dict_row = {'Date': date, 'Sensor':landsat, 'WSID':wsid, 'Upslope_NDVI':ndvi[0], 'Downslope_NDVI':ndvi[1], 'RatioDownUp_NDVI':ndvi[2], 'SD_NDVI':ndvi[3], 'Mean_NDVI':ndvi[4]}
        row_list.append(dict_row)
    
    # put everything in a dataframe
    df = pd.DataFrame(row_list)

    # save the dataframe 
    ndvi_path = os.path.join(home, "Data", "NDVI", "catchment_ratio_ndvi_results")
   # is_exist = os.path.exists(ndvi_path)
   # if not is_exist:
   #     os.makedirs(ndvi_path)
    df.to_csv(os.path.join(ndvi_path, "ratio_ndvi_" + str(wsid) + ".csv"))

    return('finished ' + str(wsid))


    # now find the trend using a simple linear regression for each catchment and put that into a dataframe with columns 
def trend_ratio_ndvi(file):
    # Calculate the annual trend in the ratio of down:up slope ndvi and the trend in SD 
    # Input: 
    #   file = csv file of a catchment with all of the ndvi info in catchment_ratio_ndvi_result
    # Output: 
    #   dictionary of the results - use dictioanry so it can be easily appended to lsit and then coereced  - supposed to be more memory / time efficient   
    if os.path.getsize(file) == 0:
        return 
    df = pd.read_csv(file)
    df['Year'] = df['Date'].astype(str).str[:4].astype(int)
    df = df[-df.eq(-9999.0).any(1)]
    df = df[-df.eq(0.0).any(1)]
    df = df.dropna()
    
    df['SD_NDVI'] = pd.to_numeric(df['SD_NDVI'])
    df['RatioDownUp_NDVI'] = pd.to_numeric(df['RatioDownUp_NDVI'])

    annual_summary = df.groupby('Year')['RatioDownUp_NDVI'].mean().to_frame().reset_index()
    slope, intercept, r_value, p_value, std_err = stats.linregress(annual_summary.Year,annual_summary.RatioDownUp_NDVI)

    annual_summary_sd = df.groupby('Year')['SD_NDVI'].mean().to_frame().reset_index()
    slope_sd, intercept_sd, r_value_sd, p_value_sd, std_err_sd = stats.linregress(annual_summary_sd.Year,annual_summary_sd.SD_NDVI)

    annual_summary_upslope = df.groupby('Year')['Upslope_NDVI'].mean().to_frame().reset_index()
    slope_up, intercept_up, r_value_up, p_value_up, std_err_up = stats.linregress(annual_summary_upslope.Year,annual_summary_upslope.Upslope_NDVI)

    annual_summary_downslope = df.groupby('Year')['Downslope_NDVI'].mean().to_frame().reset_index()
    slope_down, intercept_down, r_value_down, p_value_down, std_err_down = stats.linregress(annual_summary_downslope.Year,annual_summary_downslope.Downslope_NDVI)


    df['Mean_NDVI'] = pd.to_numeric(df['Mean_NDVI'])
    mean_ndvi = np.mean(df['Mean_NDVI'])

    wsid = df.WSID.iloc[0]
    ratio_result_dict = {'wsid':wsid, 'mean_ndvi':mean_ndvi ,'metric':'trend_ratio_ndvi', 'slope':slope, 'intercept':intercept, 'r_value':r_value, 'p_value':p_value, 'std_err':std_err, 'n':df.shape[0], 'years':annual_summary.shape[0]}
    sd_result_dict = {'wsid':wsid, 'mean_ndvi':mean_ndvi ,'metric':'trend_sd_ndvi','slope':slope_sd, 'intercept':intercept_sd, 'r_value':r_value_sd, 'p_value':p_value_sd, 'std_err':std_err_sd, 'n':df.shape[0], 'years':annual_summary.shape[0]}
    up_result_dict = {'wsid':wsid, 'mean_ndvi':mean_ndvi ,'metric':'trend_upslope_ndvi','slope':slope_up, 'intercept':intercept_up, 'r_value':r_value_up, 'p_value':p_value_up, 'std_err':std_err_up, 'n':df.shape[0], 'years':annual_summary.shape[0]}
    down_result_dict = {'wsid':wsid, 'mean_ndvi':mean_ndvi ,'metric':'trend_downslope_ndvi','slope':slope_down, 'intercept':intercept_down, 'r_value':r_value_down, 'p_value':p_value_down, 'std_err':std_err_down, 'n':df.shape[0], 'years':annual_summary.shape[0]}
    result_dict = [ratio_result_dict, sd_result_dict, up_result_dict, down_result_dict]

    return result_dict

def trend_ratio_ndvi_outliers(file):
    # Calculate the annual trend in the ratio of down:up slope ndvi and the trend in SD 
    # Input: 
    #   file = csv file of a catchment with all of the ndvi info in catchment_ratio_ndvi_result
    # Output: 
    #   dictionary of the results - use dictioanry so it can be easily appended to lsit and then coereced  - supposed to be more memory / time efficient   
    if os.path.getsize(file) == 0:
        return 
    df = pd.read_csv(file)
    df['Year'] = df['Date'].astype(str).str[:4].astype(int)
    df = df[-df.eq(-9999.0).any(1)]
    df = df[-df.eq(0.0).any(1)]
    df = df.dropna()
    
    df['SD_NDVI'] = pd.to_numeric(df['SD_NDVI'])
    df['RatioDownUp_NDVI'] = pd.to_numeric(df['RatioDownUp_NDVI'])

    # three sigma rule 
    ds_low = np.mean(df.Downslope_NDVI) - (np.std(df.Downslope_NDVI) * 3)
    ds_high = np.mean(df.Downslope_NDVI) + (np.std(df.Downslope_NDVI) * 3)
    us_low = np.mean(df.Upslope_NDVI) - (np.std(df.Upslope_NDVI) * 3)
    us_high = np.mean(df.Upslope_NDVI) + (np.std(df.Upslope_NDVI) * 3)
    ratio_low = np.mean(df.RatioDownUp_NDVI) - (np.std(df.RatioDownUp_NDVI) * 3)
    ratio_high = np.mean(df.RatioDownUp_NDVI) + (np.std(df.RatioDownUp_NDVI) * 3)
    sd_low = np.mean(df.SD_NDVI) - (np.std(df.SD_NDVI) * 3)
    sd_high = np.mean(df.SD_NDVI) + (np.std(df.SD_NDVI) * 3)  
    df = df[-((df['Upslope_NDVI'] < us_low) | (df['Upslope_NDVI'] > us_high) | (df['Downslope_NDVI'] < ds_low) | (df['Downslope_NDVI'] > ds_high) | (df['RatioDownUp_NDVI'] < ratio_low) | (df['RatioDownUp_NDVI'] > ratio_high) | (df['SD_NDVI'] < sd_low) | (df['SD_NDVI'] > sd_high))]      

    annual_summary = df.groupby('Year')['RatioDownUp_NDVI'].mean().to_frame().reset_index()
    slope, intercept, r_value, p_value, std_err = stats.linregress(annual_summary.Year,annual_summary.RatioDownUp_NDVI)

    annual_summary_sd = df.groupby('Year')['SD_NDVI'].mean().to_frame().reset_index()
    slope_sd, intercept_sd, r_value_sd, p_value_sd, std_err_sd = stats.linregress(annual_summary_sd.Year,annual_summary_sd.SD_NDVI)

    annual_summary_upslope = df.groupby('Year')['Upslope_NDVI'].mean().to_frame().reset_index()
    slope_up, intercept_up, r_value_up, p_value_up, std_err_up = stats.linregress(annual_summary_upslope.Year,annual_summary_upslope.Upslope_NDVI)

    annual_summary_downslope = df.groupby('Year')['Downslope_NDVI'].mean().to_frame().reset_index()
    slope_down, intercept_down, r_value_down, p_value_down, std_err_down = stats.linregress(annual_summary_downslope.Year,annual_summary_downslope.Downslope_NDVI)


    df['Mean_NDVI'] = pd.to_numeric(df['Mean_NDVI'])
    mean_ndvi = np.mean(df['Mean_NDVI'])

    wsid = df.WSID.iloc[0]
    ratio_result_dict = {'wsid':wsid, 'mean_ndvi':mean_ndvi ,'metric':'trend_ratio_ndvi', 'slope':slope, 'intercept':intercept, 'r_value':r_value, 'p_value':p_value, 'std_err':std_err, 'n':df.shape[0], 'years':annual_summary.shape[0]}
    sd_result_dict = {'wsid':wsid, 'mean_ndvi':mean_ndvi ,'metric':'trend_sd_ndvi','slope':slope_sd, 'intercept':intercept_sd, 'r_value':r_value_sd, 'p_value':p_value_sd, 'std_err':std_err_sd, 'n':df.shape[0], 'years':annual_summary.shape[0]}
    up_result_dict = {'wsid':wsid, 'mean_ndvi':mean_ndvi ,'metric':'trend_upslope_ndvi','slope':slope_up, 'intercept':intercept_up, 'r_value':r_value_up, 'p_value':p_value_up, 'std_err':std_err_up, 'n':df.shape[0], 'years':annual_summary.shape[0]}
    down_result_dict = {'wsid':wsid, 'mean_ndvi':mean_ndvi ,'metric':'trend_downslope_ndvi','slope':slope_down, 'intercept':intercept_down, 'r_value':r_value_down, 'p_value':p_value_down, 'std_err':std_err_down, 'n':df.shape[0], 'years':annual_summary.shape[0]}
    result_dict = [ratio_result_dict, sd_result_dict, up_result_dict, down_result_dict]
    
    return result_dict



def new_catchment_ndvi_ratio_ts(boundary, wsid, uca, ndvi_file_list):
    # Summarize upslope and downslope ndvi + the ratio for all images with no missing data and save to a csv for the given catchment
    # Inputs 
    #   boundary= catchment boundary 
    #   wsid = catchment ID 
    #   uca = usa filepath at 10m resolution for whole southern apps
    #   ndvi_file_list = list of all potential ndvi files to use 
    # Outputs 
    #   a csv is saved with the average upslope, downslope, and ratio downslope:upslope ndvi for every image with no missing data for a given catchment 

    # start a list to collect rows in form of dictionaries
    row_list = []
    K = 0
    idx = -1
    print(K)
    for image in ndvi_file_list: 
        idx +=1
        print(idx)
        # check for missing data within the boundary and if there is missing data skip to the next image
        # otherwise the function returns the image and transform 
        missing_data = check_missing_data(image, boundary['geometry'])

        if missing_data == False: 
            continue
        
        # on the first image with no missing data, resample UCA to match the NDVI
        K = K + 1
        if K == 1:
            uca_out = os.path.join(home, "Data", "Topography", "UCA", "uca" + str(wsid) + ".tif")
            resample_to_match_io(uca, uca_out, image, missing_data[0], missing_data[1], "bilinear")
            # classify upslope vs downslope pixels in the resampled UCA 
            upslope_mask, downslope_mask = classify_uca(uca_out, missing_data[0])
        
        # take the average ndvi for upslope and downslope and find the ratio 
        ndvi = calculate_ratio_ndvi(missing_data[0], upslope_mask, downslope_mask)

        # get the date 
        date = os.path.basename(image)[12:20]
        landsat = os.path.basename(image)[0:4]
        # create a dictionary 
        dict_row = {'Date': date, 'Sensor':landsat, 'WSID':wsid, 'Upslope_NDVI':ndvi[0], 'Downslope_NDVI':ndvi[1], 'RatioDownUp_NDVI':ndvi[2], 'SD_NDVI':ndvi[3]}
        row_list.append(dict_row)

    
    # put everything in a dataframe
    df = pd.DataFrame(row_list)

    # save the dataframe 
    ndvi_path = os.path.join(home, "Data", "NDVI", "catchment_ratio_ndvi_results")
    df.to_csv(os.path.join(ndvi_path, "ratio_ndvi_" + str(wsid) + ".csv"))
        
    return




def catchment_et_metrics(boundary, wsid, uca, et_file_list):
    # Summarize upslope and downslope ndvi + the ratio for all images with no missing data and save to a csv for the given catchment
    # Inputs 
    #   boundary= catchment boundary 
    #   wsid = catchment ID 
    #   dem = DEM in the original resolution cropped + buffer around the catchment 
    #   ndvi_file_list = list of all potential ET to use  
    # Outputs 
    #   a csv is saved with the average, sd, average upslope, downslope, ratio downslope:upslope ndvi for every image with no missing data for a given catchment 
   
    # create uca folder with the wsid so that we can do this for multiple watersheds at a time 
    uca_path = os.path.join(home, "Data", "Topography", "UCA","UCA_ET" + str(wsid))
    is_exist = os.path.exists(uca_path)
    if not is_exist:
        os.makedirs(uca_path)
    # crop the UCA to boundary and save it to the watershed specific uca folder 
    uca_path_name = os.path.join(uca_path, "uca.tif")
    crop_to_boundary(uca, boundary, uca_path_name)

    # start a list to collect rows in form of dictionaries
    row_list = []

    for image in et_file_list: 

        # check for missing data within the boundary and if there is missing data skip to the next image
        # otherwise the function returns the image and transform 
        missing_data = check_missing_data(image, boundary['geometry'])
        if missing_data == False:
            #print(K)
            continue

        # if there is no missing data, check if there is a UCA raster that matches the new cropped NDVI 
        # if not, resample the UCA to match the cropped NDVI 
        uca_files = glob(os.path.join(uca_path, "*.tif"))
        for file in uca_files:
            check = check_metadata(file, image, missing_data[0], missing_data[1])
            if check == True: 
                matching_uca = file
                break
            else:
                matching_uca = False
        if matching_uca == False:
            uca_out = os.path.join(uca_path, "uca" + str(len(uca_files) + 1) + ".tif")
            resample_to_match_io(os.path.join(uca_path, "uca.tif"), uca_out, image, missing_data[0], missing_data[1], "bilinear")
            matching_uca = uca_out
           
        # classify upslope vs downslope pixels in the resampled UCA 
        upslope_mask, downslope_mask = classify_uca(matching_uca, missing_data[0])
        
        # take the average ndvi for upslope and downslope and find the ratio 
        ET = calculate_ratio_et(missing_data[0], upslope_mask, downslope_mask)

        # get the date 
        date = os.path.basename(image)[52:65]
        # create a dictionary 
        dict_row = {'Date': date, 'Sensor':"ECOSTRESS", 'WSID':wsid, 'Upslope_ET':ET[0], 'Downslope_ET':ET[1], 'RatioDownUp_ET':ET[2], 'SD_ET':ET[3], 'Avg_ET':ET[4]}
        row_list.append(dict_row)
    
    # put everything in a dataframe
    df = pd.DataFrame(row_list)

    # save the dataframe 
    et_path = os.path.join(home, "Data", "Evapotranspiration", "catchment_et_metrics")
    df.to_csv(os.path.join(et_path, "ratio_et_" + str(wsid) + ".csv"))

    return('finished ' + str(wsid))


def summarize_et_catchment(filepath):
    # Calculate the average daytime and nightime ET, ratioDownslopeUpslope, and SD 
    # Input
    # filepath = csv corresponding to one hw catchment with et summarized for all clear sky images 
    # Output = a dictionary with summarized metrics. To be combined with other catchment ids later 
    if os.path.getsize(filepath)==0:
        return
    df = pd.read_csv(filepath)
    if df.shape[0] == 0:
        return
    df['Hour'] = df['Date'].apply(str).str[7:9].apply(int) # get the hour 
    df = df[(df['Hour']>=10) & (df['Hour'] <= 18)] # subset to peak daylight during the summer of 10 - 6 
    if df.shape[0] == 0: 
        return
    ratio_avg = np.mean(df.RatioDownUp_ET)
    sd_avg = np.mean(df.SD_ET)
    avg_avg = np.mean(df.Avg_ET)
    result = {'wsid':df.WSID.iloc[0], 'avg_et':avg_avg, 'sd_et':sd_avg, 'ratio_et':ratio_avg}
    return result
