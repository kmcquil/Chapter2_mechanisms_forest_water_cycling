########################################################################################################################
########################################################################################################################
## Spatial regression of the trend lateral connectivity against climatology, topography, watershed characteristics 
## 
########################################################################################################################
########################################################################################################################
library(data.table)
library(sf) # 
library(nortest) 
library(caret)
library(spdep) #
library(ape)
library(MASS)
library(gstat) # 
library(nlme)

home <- "/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling/"

########################################################################################################################
########################################################################################################################
## Create the data frame of potential predictors and response variable with one for each catchment
## Center and standardize all predictors (not response)
## Check normality of response variable 
## Check multicollinearity of the response variable 

########################################################################################################################
# Climate predictors climatology summarized for headwater catchments 
get_climatology <- function(VAR){
  dt <- fread(paste0(home, "Data/Climate/Summary/", VAR,".csv"))
  climatology <- dt[Type == "climatology", c("NHDPlusID", "Value")]
  colnames(climatology) <- c("NHDPlusID", VAR)
  return(climatology)
}

# need to add in tmin after it finishes 
tmin <- get_climatology("tmin")
tmax <- get_climatology("tmax")
prcp <- get_climatology("prcp")
vp <- get_climatology("vp")
dryday <- get_climatology("drydays")
clim_dt <- Reduce(function(x, y) merge(x, y, all=TRUE), list(tmin, tmax, prcp, vp, dryday))

########################################################################################################################
# Topography summarized for headwater catchments 
get_topo <- function(VAR){
  dt <- fread(paste0(home, "Data/Topography/", VAR,"_hw_summary.csv"))
  colnames(dt) <- c("V1", "NHDPlusID", VAR)
  return(dt <- dt[,c(2,3)])
}
topo_dt <- Reduce(function(x,y) merge(x,y, all=TRUE), list(get_topo("Aspect"),get_topo("Slope"), get_topo("Elevation"), get_topo("Latitude"), get_topo("HLI")))


########################################################################################################################
# NLCD 
# 41 = Deciduous  42 = Evergreen   43 = Mixed 
# calculate the % of the catchment that is deciduous or mixed or something 
nlcd_2001 <- fread(paste0(home, "Data/nlcd_permanent_forest/hw_nlcd_2001_summary.csv"))
nlcd_2016 <- fread(paste0(home, "Data/nlcd_permanent_forest/hw_nlcd_2016_summary.csv"))

nlcd_2001[,3:ncol(nlcd_2001)] <- nlcd_2001[,3:ncol(nlcd_2001)]/rowSums(nlcd_2001[,3:ncol(nlcd_2001)], na.rm=TRUE)
nlcd_2016[,3:ncol(nlcd_2016)] <- nlcd_2016[,3:ncol(nlcd_2016)]/rowSums(nlcd_2016[,3:ncol(nlcd_2016)], na.rm=TRUE)

# match the columns of 2001 to 2016 
setcolorder(nlcd_2016, names(nlcd_2001))
nlcd_2001[is.na(nlcd_2001)] = 0
nlcd_2016[is.na(nlcd_2016)] = 0

# find the difference between 41, 42, 43
diff <- nlcd_2016[,3:ncol(nlcd_2016)] - nlcd_2001[,3:ncol(nlcd_2001)]
diff$total <- rowSums(diff)

# % of headwater catchments that did not change  = 89% 
nrow(diff[total==0])/nrow(diff)

# % of headwater catchments that 41, 42, 43 did not change  = 87% 
diff$total_forest <- rowSums(diff[,c("41", "42", "43")])
nrow(diff[total_forest == 0])/nrow(diff)
hist(diff[!total_forest == 0]$total_forest)

# % evergreen   
nlcd_summary <- (nlcd_2001[,3:ncol(nlcd_2001)] + nlcd_2016[,3:ncol(nlcd_2016)])/2
nlcd_summary$pct_evergreen <- nlcd_summary$`42` # this code obviously isn't doing that
nlcd_summary <- cbind(nlcd_2016[,c("NHDPlusID")], nlcd_summary)

########################################################################################################################
# watershed area in km2
roi <- st_read(paste0(home, "/Data/Catchments/Headwater/headwater_catchments_perm_forest_32617.shp"))
area <- as.data.frame(roi)[,c('NHDPlusID', 'AreaSqKm_x')]

########################################################################################################################
# create a final dt of all predictor vars 
# Only use the ones with a significant trend 
pred <- Reduce(function(x,y) merge(x,y,left=TRUE), list(clim_dt, topo_dt, nlcd_summary[,c("NHDPlusID", "pct_evergreen")], area))

pred_long <- as.data.table(gather(pred, var, value, tmin:AreaSqKm_x))
summary_preds <- pred_long[,.(first_quartile = round(quantile(value, 0.25), 2), 
             median_quartile = round(quantile(value, 0.5), 2), 
             mean_quartile = round(mean(value), 2), 
             third_quartile = round(quantile(value, 0.75), 2), 
             sd_quartile = round(sd(value),2)), .(var)]
colnames(summary_preds) <- c("Predictor", "First Quartile", "Median", "Mean", "Third Quartile", "SD")
fwrite(summary_preds, paste0(home, "/Data/summary_predictors.csv"))

pred_sd <- cbind(pred[,c("NHDPlusID")], scale(pred[,2:ncol(pred)], center = TRUE, scale = TRUE))

########################################################################################################################
# Response variable joined with the spatial data 
resp <- fread(paste0(home, "Data/NDVI/ratio_ndvi_trend_results.csv"))
resp$wsid <- as.numeric(resp$wsid)
resp_sig <- resp[metric == "trend_ratio_ndvi" & p_value <= 0.05,] # keep only the responses that were significant 
resp_sig_sub <- resp_sig[,c("wsid","slope", "mean_ndvi")]
colnames(resp_sig_sub) <- c("NHDPlusID", "Trend", "Mean_NDVI")
roi_sig <- roi[roi$NHDPlusID %in% resp_sig$wsid, c("NHDPlusID")]
roi_sig <- merge(roi_sig, resp_sig_sub,by="NHDPlusID")

# final df of the response and predictors 
final_dt <- merge(roi_sig, pred_sd)

# check normality of the response variable 
hist(final_dt$Trend, breaks = 50)
ad.test(final_dt$Trend)
# according to this test the data is not normally distributed 
# however, I have over 10,000 rows so I think it's fine 

# check the multicollinearity of predictor variables
# Want VIF <= 5
lin_mod <- lm(Trend ~ Aspect + Elevation + Slope + pct_evergreen + AreaSqKm_x + Latitude, data = final_dt)
lin_mod <- lm(Trend ~ Elevation + HLI + pct_evergreen + AreaSqKm_x, data = final_dt)
car::vif(lin_mod)

cor.test(final_dt$Trend, final_dt$Elevation)
cor.test(final_dt$Trend, final_dt$Aspect)
cor.test(final_dt$Trend, final_dt$Slope)
cor.test(final_dt$Trend, final_dt$pct_evergreen)
cor.test(final_dt$Trend, final_dt$AreaSqKm_x)
cor.test(final_dt$Trend, final_dt$Latitude)


# fuck 
lin_mod <- lm(Trend ~ Aspect + Elevation + Slope + pct_evergreen + AreaSqKm_x + Latitude, data = final_dt)
lin_mod <- lm(Trend ~ Elevation + HLI + pct_evergreen + AreaSqKm_x, data = final_dt)

lin_mod <- lm(Mean_NDVI ~ Aspect + Elevation + Slope + pct_evergreen + AreaSqKm_x + Latitude, data = final_dt)
lin_mod <- lm(Mean_NDVI ~ Aspect + Elevation + Slope + pct_evergreen + AreaSqKm_x + Latitude + tmin + tmax + prcp + vp + drydays, data = final_dt)

lin_mod <- lm(Mean_NDVI ~ Elevation + HLI + pct_evergreen + AreaSqKm_x, data = final_dt)
lin_mod <- lm(Mean_NDVI ~ Elevation + HLI + pct_evergreen + AreaSqKm_x + tmin + tmax + prcp + vp + drydays, data = final_dt)


summary(lin_mod)

car::vif(lin_mod)


########################################################################################################################
########################################################################################################################
## Test for spatial autocorrelation using Moran's i 
# To calculate Moran’s I, generate a matrix of inverse distance weights
# Since the polygons are close together, use a reprojected crs (utm17N)
# find the centroids of each polygon 

#final_dt <- final_dt[1:1000,] # just for testing purposes on a smaller subset 
final_dt <- cbind(final_dt, st_coordinates(st_centroid(final_dt)))
final_dt_notshp <- st_drop_geometry(final_dt)
fwrite(final_dt_notshp, paste0(home, "/Data/SpatialRegression/final_dt.csv"))

dists <- as.matrix(dist(cbind(final_dt_notshp$X, final_dt_notshp$Y)))
inv_dists <- 1/dists
diag(inv_dists) <- 0
# We have created a matrix where each off-diagonal entry [i, j] in the matrix is equal to 1/(distance between point i and point j).
# null hypothesis = no spatial autocorrelation. so if p <= 0.05, spatial autocorrelation is present
morans <- Moran.I(final_dt_notshp$Trend, inv_dists)

########################################################################################################################
########################################################################################################################
## Fit the spatial autocorrelation 
## https://towardsdatascience.com/spatial-regression-using-fabricated-data-bbdb35da4851
# start by making a variogram 

# gstat package 
vgm <- variogram(Trend ~ Y + X, data = final_dt, cutoff = 250000, width = 100) # cutoff = 5000, width = 50
plot(vgm)

fit_vgm = fit.variogram(vgm, 
                        vgm(c("Gau", "Sph", "Mat", "Exp")), 
                        fit.kappa = TRUE, 
                        fit.method= 7)
plot(fit_vgm, cutoff = 15000)


########################################################################################################################
########################################################################################################################
# Use a generalized linear model with backward step-wise regression to find the best model 
# Do this for each spatial correlation structure and compare the AIC of each final option 

# start the hpc part of the code 
library(data.table)
library(caret)
library(MASS)
library(nlme)
library(doParallel)  #Foreach Parallel Adaptor 
library(foreach)     #Provides foreach looping construct
UseCores <- 5
clust <- makeCluster(UseCores)
registerDoParallel(clust)

home <- "/share/klmarti3/kmcquil/Chapter2_mechanisms_forest_water_cycling/"
final_dt_notshp <- fread(paste0(home, "/Data/SpatialRegression/final_dt.csv"))
spat_struc <- list(corSpher(form =~ X + Y, nugget = TRUE), corExp(form =~ X + Y, nugget = TRUE), corGaus(form =~ X + Y, nugget = TRUE), corRatio(form =~ X + Y, nugget = TRUE))

potential_models <- foreach(i = 1:length(spat_struc))%dopar%{
  
  library(nlme)
  library(MASS)
  library(caret)
  
  .GlobalEnv$final_dt_notshp <- final_dt_notshp
  .GlobalEnv$spat_cor <- spat_struc[[i]]

  model <- gls(Trend ~ Elevation + HLI + pct_evergreen + AreaSqKm_x + X + Y, 
                        data=final_dt_notshp, 
                        correlation = spat_cor,
                        method="ML")  
  
  eval_models <- MASS::stepAIC(model, direction='backward', trace=TRUE, data = final_dt_notshp)
  
}

saveRDS(potential_models, file = paste0(home, "/Data/SpatialRegression/potential_models.rds"))

potential_models <- readRDS(paste0(home, "/Data/SpatialRegression/potential_models.rds"))

# lower indicates better model fit 
AIC(potential_models[[1]])
AIC(potential_models[[2]])
AIC(potential_models[[3]])
AIC(potential_models[[4]])


plot(Variogram(potential_models[[1]]), 
     main = "Spherical Correlation")
plot(Variogram(potential_models[[2]]), 
     main = "Exponential Correlation")
plot(Variogram(potential_models[[3]]), 
     main = "Gaussian Correlation")
plot(Variogram(potential_models[[4]]), 
     main = "Quadratic Correlation")





###########################################################################################
## It's weird so little variability is explained by these predictors 
## It doesn't even seem like there is a spatial trend in the trends in the ratio 
library(corrplot)
library(RColorBrewer)
library(ggplot2)
library(ggcorrplot)

## Look at the correlation of regular NDVI 
sub_mean_ndvi <- st_drop_geometry(final_dt[,c("Trend","Mean_NDVI", "tmin", "tmax", "prcp", "vp", "drydays", "Aspect", "Elevation", "Slope", "Latitude", "HLI", "pct_evergreen", "AreaSqKm_x")])
M <-cor(sub_mean_ndvi)
ggcorrplot(M, 
           hc.order = TRUE, 
           type = "lower",
           lab = TRUE)







########################################################################################################
## updated 

final_dt_notshp[,c("X", "Y")] <- scale(final_dt_notshp[,c("X", "Y")], center=TRUE, scale=TRUE)

model_trend <- lm(Trend ~ Elevation + Slope + HLI + pct_evergreen + AreaSqKm_x + X + Y, data=final_dt_notshp)
step_model_trend <- MASS::stepAIC(model_trend, direction='backward', trace=TRUE, data = final_dt_notshp)
car::vif(step_model_trend)
# 0.08578




model_ndvi <- lm(Mean_NDVI ~ Elevation + Slope + HLI + pct_evergreen + AreaSqKm_x + Latitude + X + Y, data=final_dt_notshp)
step_model_ndvi <- MASS::stepAIC(model_ndvi, direction='backward', trace=TRUE, data = final_dt_notshp)
car::vif(step_model_ndvi)

model_ndvi <- lm(Mean_NDVI ~ Elevation + Slope + HLI + pct_evergreen + AreaSqKm_x + X + Y, data=final_dt_notshp)
step_model_ndvi <- MASS::stepAIC(model_ndvi, direction='backward', trace=TRUE, data = final_dt_notshp)
car::vif(step_model_ndvi)





##############################################################################
# group by significance and plot histograms 
library(ggplot2)
resp <- fread(paste0(home, "/Data/NDVI/ratio_ndvi_trend_results.csv"))
resp <- resp[resp$metric == "trend_ratio_ndvi",]
resp$wsid <- as.numeric(resp$wsid)
resp_sub <- resp[,c("wsid", "slope", "mean_ndvi", "p_value")]
colnames(resp_sub) <- c("NHDPlusID", "Trend", "Mean_NDVI", "Pvalue")
roi <- st_read(paste0(home, "/Data/Catchments/Headwater/headwater_catchments_perm_forest_32617.shp"))
roi_sig <- roi[roi$NHDPlusID %in% resp_sub$NHDPlusID, c("NHDPlusID")]
roi_sig <- merge(roi_sig, resp_sub,by="NHDPlusID")
final_dt <- merge(roi_sig, pred)
final_dt <- st_drop_geometry(final_dt)
fwrite(final_dt, paste0(home, "/Data/catchment_avg_vars.csv"))


