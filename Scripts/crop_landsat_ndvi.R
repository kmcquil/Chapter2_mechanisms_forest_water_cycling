library(raster)
library(doParallel)
library(foreach)


UseCores <- 5
cl <- makeCluster(UseCores)
registerDoParallel(cl)

home <- "/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling/Data/NDVI/Landsat"
ndvi_files <- list.files(paste0(home, "/Landsat_images_exported_gee"), full.names=TRUE, pattern = "*.tif$")
ndvi_names <- list.files(paste0(home, "/Landsat_images_exported_gee"), pattern = "*.tif$")

foreach(i = 1:length(ndvi_files))%dopar%{
  library(raster)
  r1 <- raster(ndvi_files[i])
  r1NaM <- is.na(as.matrix(r1))
  colNotNA <- which(colSums(r1NaM) != nrow(r1))
  rowNotNA <- which(rowSums(r1NaM) != ncol(r1))
  r3Extent <- extent(r1, rowNotNA[1], rowNotNA[length(rowNotNA)],
                     colNotNA[1], colNotNA[length(colNotNA)])
  r3 <- crop(r1, r3Extent)
  
  writeRaster(r3, paste0(home, "/Landsat_images_cropped/", ndvi_names[i]))
}

