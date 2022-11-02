########################################################################################
########################################################################################
# Calculate the heat loading index following equations from McCune 2007 
# Code based on the hli function from the SpatialEco package 
library(raster)
library(rgdal)
library(gdalUtils)
library(doParallel)  #Foreach Parallel Adaptor 
library(foreach)     #Provides foreach looping construct

#home <- "/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling/"
home <- "/share/klmarti3/kmcquil/Chapter2_mechanisms_forest_water_cycling/"

UseCores <- 20
clust <- makeCluster(UseCores)
registerDoParallel(clust)

warpMn <- function(fin, template, fout,re){
  
  res <- res(template)
  t1 <- c(xmin(template), ymin(template), 
          xmax(template), ymax(template))
  res_out <- crs(template)
  
  gdalwarp(fin, fout, t_srs=res_out, tr=res, te=t1, r=re, overwrite = T)
}
fin <- paste0(home, "Data/Topography/Elevation/elevation_10m_sbr.tif")
fout <- paste0(home, "Data/Topography/Elevation/elevation_10m_sbr_32617.tif")
gdalwarp(fin, fout, re='bilinear', t_srs="+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs", overwrite=T)
elevation <- raster(fout)

lat_fin <- paste0(home, "Data/latitude.tif")
lat_fout <- paste0(home, "Data/latitude_32617.tif")
re='bilinear'
warpMn(lat_fin, elevation, fout, re)
latitude <- raster(lat_fout)

l <- raster::clusterR(latitude, fun = function(x){abs(x)*0.017453293}, cl = clust)
cl <- clusterR(l, fun = cos, cl = clust)
sl <- clusterR(l, fun = sin, cl = clust)
#l <- raster::calc(latitude, fun=function(x){abs(x)*0.017453293})
#cl <- raster::calc(l, fun=cos)
#sl <- raster::calc(l, fun=sin)

tmp1 <- raster::clusterR(elevation, fun=raster::terrain, args=list(opt="slope", unit="degrees"), cl = clust)
tmp1 <- raster::clusterR(tmp1, fun=function(x){x*0.017453293}, cl = clust)
tmp2 <- raster::clusterR(elevation, fun=raster::terrain, args=list(opt="aspect", unit="degrees"), cl = clust)
tmp2 <- raster::clusterR(tmp2, fun=function(x){x*0.017453293}, cl = clust)
#tmp1 <- raster::terrain(elevation, opt="slope", unit="degrees") * 0.017453293
#tmp2 <- raster::terrain(elevation, opt="aspect", unit="degrees") * 0.017453293

# Folded Aspect Northern Hemisphere  (180 - (Aspect – 225) )
# 180(deg)=3.141593(rad), 225=3.92699(rad)
tmp3 <- raster::clusterR(tmp2, fun=function(x) { abs(3.141593 - abs(x - 3.926991)) }, cl = clust)
#tmp3 <- raster::calc(tmp2, fun=function(x) { abs(3.141593 - abs(x - 3.926991)) } ) 

tmp4 <- raster::clusterR(tmp1, fun=cos, cl = clust)
tmp5 <- raster::clusterR(tmp1, fun=sin, cl = clust)
tmp6 <- raster::clusterR(tmp3, fun=cos, cl = clust)
tmp7 <- raster::clusterR(tmp3, fun=sin, cl=clust)
#tmp4 <- raster::calc(tmp1, fun = cos)
#tmp5 <- raster::calc(tmp1, fun = sin)
#tmp6 <- raster::calc(tmp3, fun = cos)
#tmp7 <- raster::calc(tmp3, fun = sin)

h <- exp( -1.467 +  1.582 * cl * tmp4  - 1.5 * tmp6 * tmp5 * sl - 0.262 * sl * tmp5  + 0.607 * tmp7 * tmp5)

if(cellStats(h,"max") > 1){
  hh <- ( h / cellStats(h, "max", asSample=FALSE) ) 
}	

writeRaster(h, paste0(home, "Data/Topography/HLI.tif"), overwrite=TRUE)
writeRaster(hh, paste0(home, "Data/Topography/HLI_max.tif"), overwrite=TRUE)
