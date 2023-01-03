######################################################################################
## Phenology from Taehee 

library(R.matlab)
library(data.table)

centroids <- st_read("/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling/Data/Catchments/Headwater/headwater_catchments_perm_forest_32617_centroids.shp")
phen <- readMat('/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling/Data/GIMMS_Phenology/phenology_gimms_Katie.mat')

# get the NHDPlusID and row number lists and turn into a dt 
ids <- as.data.table(phen$NHDplusID)
colnames(ids) <- c("NHDPlusID", "PixelID")
ids$NHDPlusID <- centroids$NHDPlusID
ids$PixelID <- as.character(ids$PixelID)

# create a dt of pixel id, year, GSL, leaf off day, leaf on day 
mat_to_dt <- function(mat_obj, var_name){
  
  years <- as.data.table(t(phen$yr))
  colnames(years) <- "Year"
  
  dt <- as.data.table(mat_obj)
  colnames(dt) <- as.character(seq(1, 561))
  dt <- cbind(years, dt)
  dt_long <- melt(dt, id.vars="Year", measure.vars=colnames(dt)[2:ncol(dt)], variable.name="PixelID", value.name=var_name, variable.factor=FALSE)
  return(dt_long)
}

gsl <- mat_to_dt(phen$GSL2, "GSL")
leaf_on <- mat_to_dt(phen$mid.leafon2.yr, "LeafOn")
leaf_off <- mat_to_dt(phen$mid.leafoff2.yr, "LeafOff")
phen_metrics <- Reduce(merge, list(gsl, leaf_on, leaf_off))

# merge that back to the nhdplusid
phen_dt <- merge(phen_metrics, ids, by='PixelID', allow.cartesian = TRUE)

# 1. Calculate the trend in gsl, leaf on, leaf off in each pixel 
# 2. Calculate trend in gsl, leaf on, and leaf off in each catchment (just join the trend from the pixel and get a number)
# 3. Calculate the linear model for each catchment of ndvi ratio ~ gsl , ndvi ratio ~ leaf on, ndvi ratio ~ leaf off 
# 4. Aggregate catchments by pixel, average the ndvi, and look for a trend that way 

# 1. Calculate the trend in gsl, leaf on, leaf off in each pixel 
pix_ids <- unique(phen_metrics$PixelID)
K = 0
for(id in pix_ids){
  K = K + 1
  dt <- phen_metrics[PixelID == id,]
  gsl_trend <- lm(dt$GSL ~ dt$Year)
  leaf_on_trend <- lm(dt$LeafOn ~ dt$Year)
  leaf_off_trend <- lm(dt$LeafOff ~ dt$Year)
  
  if(K == 1){
    pixel_results <- data.table('PixelID'=id, 
                                'gsl_slope'=gsl_trend$coefficients[2], 
                                'gsl_pvalue'=summary(gsl_trend)$coefficients[2,4],
                                'leafon_slope'=leaf_on_trend$coefficients[2], 
                                'leafon_pvalue'=summary(leaf_on_trend)$coefficients[2,4],
                                'leafoff_slope'=leaf_off_trend$coefficients[2], 
                                'leafoff_pvalue'=summary(leaf_off_trend)$coefficients[2,4])
  }else{
    row <- data.table('PixelID'=id, 
                      'gsl_slope'=gsl_trend$coefficients[2], 
                      'gsl_pvalue'=summary(gsl_trend)$coefficients[2,4],
                      'leafon_slope'=leaf_on_trend$coefficients[2], 
                      'leafon_pvalue'=summary(leaf_on_trend)$coefficients[2,4],
                      'leafoff_slope'=leaf_off_trend$coefficients[2], 
                      'leafoff_pvalue'=summary(leaf_off_trend)$coefficients[2,4])
    pixel_results <- rbind(results, row)
  }
}

# get fraction of pixels with significant changes 
count_gsl <- nrow(pixel_results[gsl_slope >0 & gsl_pvalue <= 0.05,])/nrow(pixel_results)
count_leafon <- nrow(pixel_results[leafon_slope < 0 & leafon_pvalue <= 0.05,])/nrow(pixel_results)
count_leafoff <- nrow(pixel_results[leafoff_slope > 0 & leafoff_pvalue <= 0.05,])/nrow(pixel_results)

# calculate how many catchments with a significant trend in lateral connectivity are in a pixel with sig trend in phen metric
ids <- as.data.table(phen$NHDplusID)
colnames(ids) <- c("NHDPlusID", "PixelID")
ids$NHDPlusID <- centroids$NHDPlusID
ids$PixelID <- as.character(ids$PixelID)

ndvi_trend <- merge(ids, resp, by.x='NHDPlusID', by.y='wsid') # resp was just sig trends in ndvi ratio
ndvi_trend <- merge(ndvi_trend, pixel_results, by='PixelID', all.x=TRUE)
ndvi_trend <- unique(ndvi_trend)

nrow(ndvi_trend[slope < 0 & gsl_slope > 0 & gsl_pvalue <= 0.05,])/nrow(ndvi_trend) # 57%
nrow(ndvi_trend[slope < 0 & leafon_slope < 0 & leafon_pvalue <= 0.05,])/nrow(ndvi_trend) # 26% 
nrow(ndvi_trend[slope < 0 & leafoff_slope > 0 & leafoff_pvalue <= 0.05,])/nrow(ndvi_trend) # 66% 

# merge the row and column into the pixel ID 
rowcol <- as.data.table(phen$rowcol)
colnames(rowcol) <- c("row", "col")
rowcol$PixelID <- seq(1, nrow(rowcol))
rowcol$PixelID <- as.character(rowcol$PixelID)
pixel_results <- merge(pixel_results, rowcol, by="PixelID", all.x=T)
pixel_results <- unique(pixel_results)

# 2. Calculate trend in gsl, leaf on, and leaf off in each catchment (just join the trend from the pixel and get a number)
ids <- as.data.table(phen$NHDplusID)
colnames(ids) <- c("NHDPlusID", "PixelID")
ids$NHDPlusID <- centroids$NHDPlusID
ids$PixelID <- as.character(ids$PixelID)

catchment_results <- merge(ids, pixel_results, by='PixelID')
catchment_results <- unique(catchment_results)
count_gsl_catchment <- nrow(catchment_results[gsl_slope >0 & gsl_pvalue <= 0.05,])/nrow(catchment_results)
count_leafon_catchment <- nrow(catchment_results[leafon_slope < 0 & leafon_pvalue <= 0.05,])/nrow(catchment_results)
count_leafoff_catchment <- nrow(catchment_results[leafoff_slope >0 & leafoff_pvalue <= 0.05,])/nrow(catchment_results)


# 3. Calculate the linear model for each catchment of ndvi ratio ~ gsl , ndvi ratio ~ leaf on, ndvi ratio ~ leaf off 
phen_dt$Year <- as.integer(phen_dt$Year)
phen_dt <- merge(phen_dt, ndvi_dt[,..ndvi_cols_keep], by=c('NHDPlusID', 'Year'))

ids <- unique(phen_dt$NHDPlusID)
phen_vars <- c("GSL", "LeafOn", "LeafOff")
K <- 0
for(var in phen_vars){
  for(id in ids){
    K <- K+1
    dt <- phen_dt[NHDPlusID == id,] 
    cols <- c(var)
    dt[, (cols) := lapply(.SD, scale), .SDcols=cols]
    corr <- lm(phen_dt[NHDPlusID == id,get("RatioDownUp_NDVI_zscore")] ~ phen_dt[NHDPlusID == id,get(var)])
    if(K == 1){
      phen_ndvi <- data.table('NHDPlusID'=id, 
                              'Variable'=var,
                              'Slope'=corr$coefficients[2],
                              'Pvalue'=summary(corr)$coefficients[2,4], 
                              'Rsquared'=summary(corr)$r.squared)
    }else{
      row = data.table('NHDPlusID'=id, 
                       'Variable'=var,
                       'Slope'=corr$coefficients[2],
                       'Pvalue'=summary(corr)$coefficients[2,4], 
                       'Rsquared'=summary(corr)$r.squared)
      phen_ndvi <- rbind(phen_ndvi, row)
    }
  }
}

lookin <- phen_ndvi[Variable == 'GSL' & Pvalue <= 0.05, ]$Rsquared
lookin <- phen_ndvi[Variable == 'GSL' & Pvalue <= 0.05, ]$Slope
nrow(phen_ndvi[Variable == 'GSL' & Pvalue <= 0.05, ])/nrow(phen_ndvi[Variable == "GSL",])
boxplot(lookin)
quantile(lookin)

lookin <- phen_ndvi[Variable == 'LeafOn' & Pvalue <= 0.05, ]$Rsquared
lookin <- phen_ndvi[Variable == 'LeafOn' & Pvalue <= 0.05, ]$Slope
nrow(phen_ndvi[Variable == 'LeafOn' & Pvalue <= 0.05, ])/nrow(phen_ndvi[Variable == "LeafOn",])
boxplot(lookin)
quantile(lookin)

lookin <- phen_ndvi[Variable == 'LeafOff' & Pvalue <= 0.05, ]$Rsquared
lookin <- phen_ndvi[Variable == 'LeafOff' & Pvalue <= 0.05, ]$Slope
nrow(phen_ndvi[Variable == 'LeafOff' & Pvalue <= 0.05, ])/nrow(phen_ndvi[Variable == "LeafOff",])
boxplot(lookin)
quantile(lookin)


# 4. Aggregate catchments by pixel, average the ndvi, and look for a trend that way 
phen_dt_pix <- phen_dt[, .(GSL=mean(GSL), 
                           LeafOn=mean(LeafOn), 
                           LeafOff=mean(LeafOff), 
                           RatioDownUp_NDVI_zscore = mean(RatioDownUp_NDVI_zscore)), .(PixelID, Year)]


ids <- unique(phen_dt_pix$PixelID)
phen_vars <- c("GSL", "LeafOn", "LeafOff")
K <- 0
for(var in phen_vars){
  for(id in ids){
    
    K <- K+1
    dt <- phen_dt_pix[PixelID == id,] 
    cols <- c(var)
    dt[, (cols) := lapply(.SD, scale), .SDcols=cols]
    
    corr <- lm(dt[,get("RatioDownUp_NDVI_zscore")] ~ dt[,get(var)])
    if(K == 1){
      phen_ndvi_pix <- data.table('PixelID'=id, 
                                  'Variable'=var,
                                  'Slope'=corr$coefficients[2],
                                  'Pvalue'=summary(corr)$coefficients[2,4], 
                                  'Rsquared'=summary(corr)$r.squared)
    }else{
      row = data.table('PixelID'=id, 
                       'Variable'=var,
                       'Slope'=corr$coefficients[2],
                       'Pvalue'=summary(corr)$coefficients[2,4], 
                       'Rsquared'=summary(corr)$r.squared)
      phen_ndvi_pix <- rbind(phen_ndvi_pix, row)
    }
  }
}

lookin <- phen_ndvi_pix[Variable == 'GSL' & Pvalue <= 0.05, ]$Rsquared
lookin <- phen_ndvi_pix[Variable == 'GSL' & Pvalue <= 0.05, ]$Slope
nrow(phen_ndvi_pix[Variable == 'GSL' & Pvalue <= 0.05, ])/nrow(phen_ndvi_pix[Variable == "GSL",])
boxplot(lookin)
quantile(lookin)

lookin <- phen_ndvi_pix[Variable == 'LeafOn' & Pvalue <= 0.05, ]$Rsquared
lookin <- phen_ndvi_pix[Variable == 'LeafOn' & Pvalue <= 0.05, ]$Slope
nrow(phen_ndvi_pix[Variable == 'LeafOn' & Pvalue <= 0.05, ])/nrow(phen_ndvi_pix[Variable == "LeafOn",])
boxplot(lookin)
quantile(lookin)

lookin <- phen_ndvi_pix[Variable == 'LeafOff' & Pvalue <= 0.05, ]$Rsquared
lookin <- phen_ndvi_pix[Variable == 'LeafOff' & Pvalue <= 0.05, ]$Slope
nrow(phen_ndvi_pix[Variable == 'LeafOff' & Pvalue <= 0.05, ])/nrow(phen_ndvi_pix[Variable == "LeafOff",])
boxplot(lookin)
quantile(lookin)






# plot the catchment results (trends of gsl, leaf on, leaf off)
plot_phen_trend <- function(VAR, PVALUE, main){
  
  df <- catchment_results[get(PVALUE) <= 0.05,]
  #lower_limit <-quantile(df[,get(VAR)], 0.01)
  #upper_limit <-quantile(df[,get(VAR)], 0.99)
  #df[get(VAR) < lower_limit, get(VAR) := lower_limit]
  
  catchment_results_sf <- merge(catch, df, by="NHDPlusID", all.y=TRUE, all.x=FALSE)
  catchment_results_cent <- cbind(catchment_results_sf, st_coordinates(st_centroid(catchment_results_sf)))
  
  
  phen_plot <- ggplot() + 
    geom_sf(data=roi, fill='white', color='black')+ 
    geom_point(data=catchment_results_cent, aes(x=X, y=Y, color=get(VAR)), alpha=0.4, size=0.001) + 
    #scale_color_gradient2(low = "dark red", mid = "light yellow", high = "dark blue", midpoint = 0, name=VAR, limits=c(lower_limit, upper_limit)) + 
    scale_color_gradient2(low = "dark red", mid = "light yellow", high = "dark blue", midpoint = 0, name=VAR) + 
    theme_classic() + 
    xlab("") + 
    ylab("") + 
    ggtitle(main) +
    theme(legend.title = element_blank()) +
    theme(legend.position = c(0.9, 0.3)) 
  #theme(legend.position = "none")
  
  return(phen_plot)
}


gsl_trend <- plot_phen_trend('gsl_slope', 'gsl_pvalue', 'GSL Trend')
leafon_trend <- plot_phen_trend('leafon_slope', 'leafon_pvalue', 'Leaf On Date Trend')
leafoff_trend <- plot_phen_trend('leafoff_slope', 'leafoff_pvalue', 'Leaf Off Date Trend')
ggarrange(gsl_trend, leafon_trend, leafoff_trend, nrow=1)

# plot the catchment results of the slope and Rsquared 
plot_phen_rels <- function(VAR){
  
  df <- phen_ndvi[Variable == VAR & Pvalue<= 0.05,]
  catchment_results_sf <- merge(catch, df, by="NHDPlusID", all.y=TRUE, all.x=FALSE)
  catchment_results_cent <- cbind(catchment_results_sf, st_coordinates(st_centroid(catchment_results_sf)))
  
  
  slope_plot <- ggplot()+ 
    geom_sf(data=roi, fill='white', color='black')+ 
    geom_point(data=catchment_results_cent, aes(x=X, y=Y, color=Slope), alpha=0.4, size=0.001) + 
    #scale_color_gradient2(low = "dark red", mid = "light yellow", high = "dark blue", midpoint = 0, name=VAR, limits=c(lower_limit, upper_limit)) + 
    scale_color_gradient2(low = "dark red", mid = "light yellow", high = "dark blue", midpoint = 0, name=VAR) + 
    theme_classic() + 
    xlab("") + 
    ylab("") + 
    ggtitle(paste0("NDVIDown:Up ~ ",VAR, " Slope"))+
    theme(legend.title = element_blank()) + 
    theme(legend.position = c(0.9, 0.3)) 
  
  r_plot <- ggplot()+ 
    geom_sf(data=roi, fill='white', color='black')+ 
    geom_point(data=catchment_results_cent, aes(x=X, y=Y, color=Rsquared), alpha=0.4, size=0.001) + 
    scale_color_viridis_c() + 
    theme_classic() + 
    xlab("") + 
    ylab("") + 
    ggtitle(paste0("NDVIDown:Up ~ ",VAR, " R2"))+
    theme(legend.title = element_blank()) + 
    theme(legend.position = c(0.9, 0.3)) 
  results <- list(slope_plot, r_plot)
  return(results)
  
}

gsl_mod <- plot_phen_rels('GSL')
on_mod <- plot_phen_rels('LeafOn')
off_mod <- plot_phen_rels('LeafOff')

ggarrange(gsl_mod[[1]], gsl_mod[[2]], 
          on_mod[[1]], on_mod[[2]], 
          off_mod[[1]], off_mod[[2]], nrow = 3, ncol=2)



########################################################################################
## Resample annual growing season, mam, son, jja tmin and vp to match gimms 
## Then pull out the specific pixel ids based on row and col from taehee and calculate pixel wise correlations 
## Make a map of the Trend in growing season length, leaf on, leaf off 
## Make a map of the correlations between growing season tmin and vp, growing season tmin and gsl, tmin and leaf on, tmin and leaf off 
library(R.matlab)
library(data.table)

phen <- readMat('/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling/Data/GIMMS_Phenology/phenology_gimms_Katie.mat')
temp_rast <- raster("/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling/Data/GIMMS_Phenology/geo00apr15a_ndvi.asc")
crs(temp_rast) <- CRS("EPSG:4326")
roi <- st_read("/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling/Data/ROI/blue_ridge_plus_reference.shp")
roi <- st_transform(roi, crs="EPSG:32617")

# resample climate anomalies to match the gimms data 
#climate_anomalies <- list.files("/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling/Data/Climate/tifs", full.names=T, pattern="*.tif$")
#climate_anomalies <- climate_anomalies[grep("anomalies_tmin_SON|anomalies_vp_SON|anomalies_tmin_MAM|anomalies_vp_MAM|anomalies_tmin_AMJJA|anomalies_vp_AMJJA", climate_anomalies)]
#climate_anomalies <- climate_anomalies[-grep("anomalies_vp_MAMJJA|anomalies_tmin_MAMJJA", climate_anomalies)]

climate_anomalies <- list.files("/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling/Data/Climate/Phenology/tifs", full.names=T, pattern="*.tif$")
for(i in climate_anomalies){
  r <- raster(i)
  r <- projectRaster(r, crs=crs(temp_rast), method='bilinear')
  r <- resample(r, temp_rast, method='bilinear')
  writeRaster(r, paste0("/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling/Data/Climate/GIMMS2/", basename(i)), overwrite = T)
}

# associate the pixel id with row/col in a table 
# loop through each tif and get the - variable, season, year, and then go to the row/col and get the value 
# put all of that into a big dataframe 
gimms_clim <- list.files("/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling/Data/Climate/GIMMS2", full.names=T)
rowcol <- as.data.table(phen$rowcol)
colnames(rowcol) <- c("row", "col")
rowcol$PixelID <- seq(1, nrow(rowcol))

library(doParallel)
cl <- makeCluster(5)
registerDoParallel(cl)

gimms_climate <- foreach(file_num = 1:length(gimms_clim), .combine=rbind)%dopar%{
  library(raster)
  library(data.table)
  
  file <- gimms_clim[file_num]
  file_strs <- strsplit(basename(file), "_")[[1]]
  var <- file_strs[2]
  season <- file_strs[3]
  year <- substr(file_strs[4], 1, 4)
  
  r <- raster(file)
  for(i in 1:nrow(rowcol)){
    PixelID <- rowcol$PixelID[i]
    anom <- r[rowcol$row[i], rowcol$col[i]]
    row <- data.table(PixelID=PixelID, 
                      Var=var, 
                      Season=season,
                      Year=year,
                      Anom=anom)
    if(i == 1){
      result_clim <- row
    }else{
      result_clim <- rbind(result_clim, row)
    }
  }
  
  result_clim
}
gimms_climate$Year <- as.numeric(gimms_climate$Year)
gimms_climate$PixelID <- as.character(gimms_climate$PixelID)



# make a dt of growing season metrics 
mat_to_dt <- function(mat_obj, var_name){
  
  years <- as.data.table(t(phen$yr))
  colnames(years) <- "Year"
  
  dt <- as.data.table(mat_obj)
  colnames(dt) <- as.character(seq(1, 561))
  dt <- cbind(years, dt)
  dt_long <- melt(dt, id.vars="Year", measure.vars=colnames(dt)[2:ncol(dt)], variable.name="PixelID", value.name=var_name, variable.factor=FALSE)
  return(dt_long)
}

gsl <- mat_to_dt(phen$GSL2, "GSL")
leaf_on <- mat_to_dt(phen$mid.leafon2.yr, "LeafOn")
leaf_off <- mat_to_dt(phen$mid.leafoff2.yr, "LeafOff")
phen_metrics <- Reduce(merge, list(gsl, leaf_on, leaf_off))

phen_metric_cors <- phen_metrics[, .(R2_LeafOn = summary(lm(GSL~LeafOn))$r.squared, 
                                     R2_LeafOff = summary(lm(GSL~LeafOff))$r.squared), .(PixelID)]

gimms_dt <- merge(gimms_climate, phen_metrics, by=c("Year", "PixelID"), all.x=TRUE, allow.cartesian=TRUE)
gimms_dt <- gimms_dt[complete.cases(gimms_dt),]

## Make a map of the correlations between growing season tmin and vp, growing season tmin and gsl, tmin and leaf on, tmin and leaf off 
## Calculate the correlation between: 
### Tmin and GSL for MAM, AMJJA, SON 
### Tmin and Leaf on for MAM, AMJJA
### Tmin and Leaf off for AMJJA, SON

cor_func <- function(X, Y){
  mod <- lm(Y ~ X)
  row <- data.table(Slope = mod$coefficients[2], 
                    Pvalue=summary(mod)$coefficients[2,4],
                    Rsquared=summary(mod)$r.squared)
  return(row)
}

PixelIDs <- unique(gimms_dt$PixelID)
Seasons <- unique(gimms_dt$Season)
Metrics <- c("GSL", "LeafOn", "LeafOff")
K <- 0
for(id in PixelIDs){
  for(season in Seasons){
    for(metric in Metrics){
      K <- K + 1
      row <- data.table(PixelID=id, X='tmin', Y=metric, Season=season)
      Y <- as.numeric(scale(gimms_dt[PixelID == id & Var == 'tmin' & Season == season, Anom], scale=T, center=T))
      X <- as.numeric(scale(gimms_dt[PixelID == id & Var == 'tmin' & Season == season, get(metric)], scale=T, center=T))
      linmod <- cbind(row, cor_func(X, Y))
      if(K == 1){
        gimms_phen_dt <- linmod
      }else{
        gimms_phen_dt <- rbind(gimms_phen_dt, linmod)
      }
    }
  }
}

par(mfrow = c(2,3))
hist(gimms_phen_dt[Season=="March" & Y == "LeafOn",]$Rsquared, xlab="R-Squared", main="Leaf On ~ Mean March Min Temp", col='#9AC791')
hist(gimms_phen_dt[Season=="April" & Y == "LeafOn",]$Rsquared, xlab="R-Squared", main="Leaf On ~ Mean April Min Temp", col='#9AC791')
hist(gimms_phen_dt[Season=="May" & Y == "LeafOn",]$Rsquared, xlab="R-Squared", main="Leaf On ~ Mean May Min Temp", col='#9AC791')

hist(gimms_phen_dt[Season=="August" & Y == "LeafOff",]$Rsquared, xlab="R-Squared", main="Leaf Off ~ Mean August Min Temp", col='#9AC791')
hist(gimms_phen_dt[Season=="September" & Y == "LeafOff",]$Rsquared, xlab="R-Squared", main="Leaf Off ~ Mean September Min Temp", col='#9AC791')
hist(gimms_phen_dt[Season=="October" & Y == "LeafOff",]$Rsquared, xlab="R-Squared", main="Leaf Off ~ Mean October Min Temp", col='#9AC791')


head(gimms_phen_dt)
highest_r2 <- gimms_phen_dt[gimms_phen_dt[, .I[Rsquared==max(Rsquared)], by=.(PixelID, Y)]$V1]
par(mfrow=c(1,2))
hist(highest_r2[Y=='LeafOn',]$Rsquared, xlab='R2', main="Leaf On ~ Mean Min Temp", ylab='# Pixels', col='#9AC791')
hist(highest_r2[Y=='LeafOff',]$Rsquared, xlab='R2', main="Leaf Off ~ Mean Min Temp",ylab='# Pixels', col='#9AC791')


# Plot a map of the trend in GSL, Leaf On, and Leaf off with black dots for significance 
# Boxplots with the R2 of pixel wise models (Tmin ~ VP), (Leaf On ~ Tmin), (Leaf off ~ Tmin), (GSL ~ Leaf on), (GSL~Leaf off)
# Bar plot of the fraction of catchments located in a pixel with sig increasing gsl, earlier leaf on, or later leaf off

# Plot the map of the trends in GSL, Leaf On, and Leaf off iwth black dots for significance 
# Get the template raster, make everything NA, and then populate the raster by looping 
temp_rast <- raster("/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling/Data/GIMMS_Phenology/geo00apr15a_ndvi.asc")
crs(temp_rast) <- CRS("EPSG:4326")
roi <- st_read("/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling/Data/ROI/blue_ridge_plus_reference.shp")
roi <- st_transform(roi, crs="EPSG:4326")

temp_rast[temp_rast < 0 | temp_rast > 0 | temp_rast == 0] <- NA
gsl_trend_rast <- temp_rast
leafon_trend_rast <- temp_rast
leafoff_trend_rast <- temp_rast

cells <- c()
for(i in 1:nrow(pixel_results)){
  cells <- c(cells, cellFromRowCol(temp_rast, pixel_results$row[i], pixel_results$col[i]))
  gsl_trend_rast[pixel_results$row[i], pixel_results$col[i]] <- pixel_results$gsl_slope[i]
  leafon_trend_rast[pixel_results$row[i], pixel_results$col[i]] <- pixel_results$leafon_slope[i]
  leafoff_trend_rast[pixel_results$row[i], pixel_results$col[i]] <- pixel_results$leafoff_slope[i]
}

gsl_trend_rast_sbr <- crop(gsl_trend_rast, roi)
leafon_trend_rast_sbr <- crop(leafon_trend_rast, roi)
leafoff_trend_rast_sbr <- crop(leafoff_trend_rast, roi)

rast_cent <- cbind(pixel_results, xyFromCell(temp_rast, cells))
rast_cent <- st_as_sf(rast_cent, coords=c('x', 'y'), crs=crs(temp_rast))

all_slopes <- c(rast_cent[rast_cent$gsl_pvalue<=0.05,]$gsl_slope, rast_cent[rast_cent$leafon_pvalue<=0.05,]$leafon_slope, rast_cent[rast_cent$leafoff_pvalue<=0.05,]$leafoff_slope)
upper_limit <- quantile(all_slopes, 0.99)
lower_limit <- quantile(all_slopes, 0.01)

test <- gsl_trend_rast_sbr
test_spdf <- as(test, "SpatialPixelsDataFrame")
test_df <- as.data.frame(test_spdf)
colnames(test_df) <- c("value", "x", "y")
test_df[test_df$value < lower_limit,]$value <- lower_limit
test_df[test_df$value > upper_limit,]$value <- upper_limit
gsl_trend_plot <- ggplot() + 
  geom_tile(data=test_df, aes(x=x, y=y, fill=value), alpha=0.8) + 
  scale_fill_gradient2(low = "red", mid = "light yellow", high = "blue", midpoint = 0,limits = c(lower_limit, upper_limit), name="") + 
  geom_sf(data=roi, fill=NA, color='black') +
  geom_sf(data=rast_cent[rast_cent$gsl_pvalue<=0.05,], color='black', size=0.1) + 
  #theme_map() +
  theme_void() + 
  theme(legend.title = element_blank()) + 
  theme(legend.position = "none") +
  #theme(legend.position="bottom") +
  #theme(legend.key.width=unit(2, "cm")) + 
  ggtitle("Growing season length") + 
  theme(plot.title = element_text(hjust = 0.5))
gsl_trend_plot

test <- leafon_trend_rast_sbr
test_spdf <- as(test, "SpatialPixelsDataFrame")
test_df <- as.data.frame(test_spdf)
colnames(test_df) <- c("value", "x", "y")
test_df[test_df$value < lower_limit,]$value <- lower_limit
test_df[test_df$value > upper_limit,]$value <- upper_limit
leafon_trend_plot <- ggplot() + 
  geom_tile(data=test_df, aes(x=x, y=y, fill=value), alpha=0.8) + 
  scale_fill_gradient2(low = "red", mid = "light yellow", high = "blue", midpoint = 0,limits = c(lower_limit, upper_limit), name="") + 
  geom_sf(data=roi, fill=NA, color='black') +
  geom_sf(data=rast_cent[rast_cent$leafon_pvalue<=0.05,], color='black', size=0.1) + 
  #theme_map() +
  theme_void() + 
  theme(legend.title = element_blank()) + 
  theme(legend.position = "none") +
  #theme(legend.position="bottom") +
  #theme(legend.key.width=unit(2, "cm")) + 
  ggtitle("Green up") + 
  theme(plot.title = element_text(hjust = 0.5))
leafon_trend_plot

test <- leafoff_trend_rast_sbr
test_spdf <- as(test, "SpatialPixelsDataFrame")
test_df <- as.data.frame(test_spdf)
colnames(test_df) <- c("value", "x", "y")
test_df[test_df$value < lower_limit,]$value <- lower_limit
test_df[test_df$value > upper_limit,]$value <- upper_limit
leafoff_trend_plot <- ggplot() + 
  geom_tile(data=test_df, aes(x=x, y=y, fill=value), alpha=0.8) + 
  scale_fill_gradient2(low = "red", mid = "light yellow", high = "blue", midpoint = 0,limits = c(lower_limit, upper_limit), name="") + 
  geom_sf(data=roi, fill=NA, color='black') +
  geom_sf(data=rast_cent[rast_cent$leafoff_pvalue<=0.05,], color='black', size=0.1) + 
  #theme_map() +
  theme_void() + 
  theme(legend.title = element_blank()) + 
  theme(legend.position = "none") +
  #theme(legend.position="bottom") +
  #theme(legend.key.width=unit(2, "cm")) + 
  ggtitle("Senesence") + 
  theme(plot.title = element_text(hjust = 0.5))
leafoff_trend_plot

legend <- ggplot() + 
  geom_tile(data=test_df, aes(x=x, y=y, fill=value), alpha=0.8) + 
  scale_fill_gradient2(low = "red", mid = "light yellow", high = "blue", midpoint = 0,limits = c(lower_limit, upper_limit), name="Days/year") + 
  geom_sf(data=roi, fill=NA, color='black') +
  geom_sf(data=rast_cent[rast_cent$leafoff_pvalue<=0.05,], color='black', size=0.1) + 
  theme_void() +
  theme(legend.position = "right", 
        legend.text=element_text(size=12)) +
  theme(axis.text=element_text(size=12, color = 'black'),
        axis.title=element_text(size=12,face="bold")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1)) + 
  theme(legend.key.width=unit(1, "cm"))
legend <- ggpubr::get_legend(legend)


ggarrange(leafon_trend_plot, leafoff_trend_plot,gsl_trend_plot, nrow=1, ncol=3, 
          common.legend = TRUE,
          legend.grob=legend,
          legend = "right")




# Boxplots with the R2 of pixel wise models (Tmin ~ VP), (Leaf On ~ Tmin), (Leaf off ~ Tmin), (GSL ~ Leaf on), (GSL~Leaf off)
# make a df with a column for pixelID, comcbo, r2

tmin_vp_r2 <- merge(tmin_vp_r2, ids, by='NHDPlusID', all.x=T)
tmin_vp_r2_pix <- tmin_vp_r2[, .(R2=mean(R2)), .(PixelID)]
tmin_vp_r2_pix$Label <- rep("VP ~ Tmin", nrow(tmin_vp_r2_pix))

highest_r2$Label <- paste0(highest_r2$Y, "~", highest_r2$X)
highest_r2_pix <- highest_r2[,c("PixelID", "Rsquared", "Label")]
colnames(highest_r2_pix) <- c("PixelID", "R2", "Label")
highest_r2_pix <- highest_r2_pix[!Label == "GSL~tmin",]
highest_r2_pix[Label == "LeafOn~tmin",Label:= "Green Up ~ Tmin"] 
highest_r2_pix[Label == "LeafOff~tmin",Label:= "Senesence ~ Tmin"]

colnames(phen_metric_cors) <- c("PixelID", "GSL ~ Green Up", "GSL ~ Senesence")
phen_metric_cors_pix <- melt(phen_metric_cors, id.vars="PixelID", measure.vars=c("GSL ~ Green Up", "GSL ~ Senesence"), 
                             variable.name="Label", value.name="R2")
phen_metric_cors_pix <- phen_metric_cors_pix[,c(1,3,2)]

#all_cors <- rbind(tmin_vp_r2_pix, highest_r2_pix, phen_metric_cors_pix)
#all_cors$Label <- factor(all_cors$Label, levels=c("Green Up ~ Tmin",  "Senesence ~ Tmin", "GSL ~ Green Up", "GSL ~ Senesence", "VP ~ Tmin"))
all_cors <- rbind(highest_r2_pix, phen_metric_cors_pix)
cors_plot <- ggplot(all_cors, aes(x=Label, y=R2, fill=Label)) + 
  geom_boxplot() + 
  #scale_fill_manual(values=c("#beaed4", "#7570b3", "#fdc086", "#bf5b17", "#a6cee3")) + 
  scale_fill_manual(values=c("#beaed4", "#7570b3", "#fdc086", "#bf5b17")) + 
  theme_classic() + 
  ylab("R-Squared") + 
  xlab("Pixel-wise Linear Model") + 
  theme(legend.position="none") + 
  theme(axis.text.x=element_text(angle=25,hjust=1))



# Bar plot of the fraction of catchments located in a pixel with sig increasing gsl, earlier leaf on, or later leaf off
ids <- as.data.table(phen$NHDplusID)
colnames(ids) <- c("NHDPlusID", "PixelID")
ids$NHDPlusID <- centroids$NHDPlusID
ids$PixelID <- as.character(ids$PixelID)

ndvi_trend <- merge(ids, resp, by.x='NHDPlusID', by.y='wsid') # resp was just sig trends in ndvi ratio
ndvi_trend <- merge(ndvi_trend, pixel_results, by='PixelID', all.x=TRUE)
ndvi_trend <- unique(ndvi_trend)

frac_gsl <- nrow(ndvi_trend[slope < 0 & gsl_slope > 0 & gsl_pvalue <= 0.05,])/nrow(ndvi_trend) # 57%
frac_leafon <- nrow(ndvi_trend[slope < 0 & leafon_slope < 0 & leafon_pvalue <= 0.05,])/nrow(ndvi_trend) # 26% 
frac_leafoff <- nrow(ndvi_trend[slope < 0 & leafoff_slope > 0 & leafoff_pvalue <= 0.05,])/nrow(ndvi_trend) # 66% 

frac_dt <- data.table(Label=c('GSL', 'Green Up', 'Senesence'), 
                      Fraction=c(round(frac_gsl*100, 2), round(frac_leafon*100, 2),round(frac_leafoff*100, 2)))
frac_dt$Label <-factor(frac_dt$Label, levels=c("Green Up", "Senesence", "GSL"))
frac_plot <- ggplot(frac_dt, aes(x=Label, y=Fraction)) + 
  geom_bar(stat='identity', aes(fill=Label)) + 
  scale_fill_manual(values=c("#7570b3", "#bf5b17", "#a6cee3")) + 
  theme_classic() + 
  theme(legend.position="none") + 
  xlab("Growing Season Metric") + 
  ylab("% of catchments") + 
  theme(axis.text.x=element_text(angle=25,hjust=1))
frac_plot


hmm <- ggarrange(leafon_trend_plot, leafoff_trend_plot,gsl_trend_plot, 
                 nrow=1, ncol=3, 
                 common.legend = TRUE,
                 legend.grob=legend,
                 legend = "right", 
                 labels='AUTO')
hmm1 <- ggarrange(cors_plot, frac_plot, nrow=1, ncol=3, labels=c("D", "E"))
hmm2 <- ggarrange(hmm, hmm1, nrow=2, ncol=1)













## test how many catchments have a significant trend in just the last 21 years 
# standardize the ratio ndvi and mean NDVI in new columns
ndvi_dt <- fread(paste0(home, "/Data/NDVI/annual_ratio_ndvi.csv"))
colnames(ndvi_dt) <- c("Year", "RatioDownUp_NDVI", "Mean_NDVI","NHDPlusID")
ndvi_dt$NHDPlusID <- as.numeric(ndvi_dt$NHDPlusID)

get_trend <- function(X, Y){
  result <- list()
  mod <- lm(Y~X)
  result[[1]] <- mod$coefficients[2]
  result[[2]] <- summary(mod)$coefficients[2,4]
  names(result) <- c("Slope", "Pvalue")
  return(result)
}

trends_1984 <- ndvi_dt[Year >=1984, get_trend(Year, RatioDownUp_NDVI), .(NHDPlusID)]
trends_2000 <- ndvi_dt[Year >=2000, get_trend(Year, RatioDownUp_NDVI), .(NHDPlusID)]
trends_2005 <- ndvi_dt[Year >=2005, get_trend(Year, RatioDownUp_NDVI), .(NHDPlusID)]

trends_pre_2000 <- ndvi_dt[Year <2000, get_trend(Year, RatioDownUp_NDVI), .(NHDPlusID)]


nrow(trends_1984[Pvalue <= 0.05 & Slope < 0,])/nrow(trends_1984)
nrow(trends_2000[Pvalue <= 0.05 & Slope < 0,])/nrow(trends_2000)
nrow(trends_2005[Pvalue <= 0.05 & Slope < 0,])/nrow(trends_2005)

nrow(trends_pre_2000[Pvalue <= 0.05 & Slope < 0,])/nrow(trends_pre_2000)










