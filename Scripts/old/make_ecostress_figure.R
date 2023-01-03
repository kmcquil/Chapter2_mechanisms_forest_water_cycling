# make a summer composite of ECOSTRESS instantaneous ET during the day 

library(raster)
library(rgdal)
library(data.table)

roi <- st_read(paste0(home, "/Data/Catchments/Reference/gages_ii/reference_keep.shp"))
ws <- roi[roi$GAGE_ID == "0344894205",]
ws <- roi[roi$GAGE_ID == "02149000",]

files <- list.files("/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling/Data/Evapotranspiration/Ecostress/Data/clean", full.names=TRUE, pattern="*.tif$")
et_files_dt <- data.table(file= files, 
                          file_short = basename(files)) 
et_files_dt$date <- substr(et_files_dt$file_short, 53, 59)
et_files_dt$hour <- substr(et_files_dt$file_short, 60, 61)

# subset to just between 12pm and 5pm 
et_day <- et_files_dt[hour >= 12 & hour <= 17,]

r <- raster(et_day$file[1])
ws <- st_transform(ws, crs=crs(r))

for(i in 1:nrow(et_day)){
  print(i)
  r <- raster(et_day$file[i])
  
  skip_to_next <- FALSE
  tryCatch(r_crop <- mask(crop(r, ws), ws), error = function(e) { skip_to_next <<- TRUE})
  print(skip_to_next)
  if(skip_to_next) { next }  
  
  r_crop <- mask(crop(r, ws), ws)
  
  if(min(values(r_crop), na.rm=T) == Inf){
    print(min(values(r_crop), na.rm=T))
    next
  }
  
  if(i == 1){
    temp <- r_crop
    writeRaster(r_crop, paste0("/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling/Data/Evapotranspiration/Ecostress/Data/clean/ws_crop/", et_day$file_short[i]), overwrite=T)
  }else{
    skip_to_next <- FALSE
    tryCatch(resample(r_crop, temp, method='bilinear'), error = function(e) { skip_to_next <<- TRUE})
    
    print(skip_to_next)
    
    if(skip_to_next) { next }  
    r_resamp <- resample(r_crop, temp, method='bilinear')
    writeRaster(r_resamp, paste0("/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling/Data/Evapotranspiration/Ecostress/Data/clean/ws_crop/", et_day$file_short[i]), overwrite=T)
  }
  
}


cropped_files <- list.files("/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling/Data/Evapotranspiration/Ecostress/Data/clean/ws_crop", full.names=T, pattern="*.tif$")
s <- raster::stack(cropped_files)
s_mean <- mean(s, na.rm=T)
plot(s_mean)
plot(st_geometry(ws),col=rgb(red = 0, green = 0, blue = 0, alpha = 0), add=T)
writeRaster(s_mean, "/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling/Data/Evapotranspiration/Ecostress/Data/clean/ws_crop/eco_mean.tif")


