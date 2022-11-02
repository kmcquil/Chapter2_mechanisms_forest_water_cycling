##################################################################################################
# Calculate a regression for each catchment x climate variable 
# Climate variables summarized as average or total for the growing season (April - August)
# This matches the NDVI summary from June - August 

library(data.table)
library(colorspace)
library(viridis)
library(RColorBrewer)
library(sf) 
library(ggplot2)
library(cowplot)
library(ggpubr)
library(patchwork)
home <- "/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling/"

#################################################################################################################################################
#################################################################################################################################################
#################################################################################################################################################
## Look at NDVI trends 

summarize_trends <- function(dt){
  frac_neg_trend <- nrow(dt[p_value <= 0.05 & slope < 0, ])/nrow(dt)
  frac_pos_trend <- nrow(dt[p_value <= 0.05 & slope > 0, ])/nrow(dt)
  frac_nonsig_trend <- 1-(frac_neg_trend + frac_pos_trend)
  return(c(frac_neg_trend, frac_pos_trend, frac_nonsig_trend))
}
ndvi_trends <- fread(paste0(home, "/Data/NDVI/ratio_ndvi_trend_results_outliers.csv"))
ratio_trends <- ndvi_trends[metric == "trend_ratio_ndvi"]
sd_trends <- ndvi_trends[metric=='trend_sd_ndvi']
upslope_trends <- ndvi_trends[metric=='trend_upslope_ndvi']
downslope_trends <- ndvi_trends[metric=='trend_downslope_ndvi']

summarize_trends(ratio_trends)
summarize_trends(sd_trends)
summarize_trends(upslope_trends)
summarize_trends(downslope_trends)

# of the catchments with a significant decreasing trend in SD or ratio NDVI, how many have a decreasing trend in downslope ndvi 
sum(downslope_trends[p_value <= 0.05 & slope < 0, ]$wsid %in% ratio_trends[p_value <= 0.05 & slope < 0, ]$wsid)/length(ratio_trends[p_value <= 0.05 & slope < 0, ]$wsid) # 11% 
sum(downslope_trends[p_value <= 0.05 & slope < 0, ]$wsid %in% sd_trends[p_value <= 0.05 & slope < 0, ]$wsid)/length(sd_trends[p_value <= 0.05 & slope < 0, ]$wsid) # 3%

# of the catchments with a significant decreasing trend in SD or ratio NDVI, how many have a decreasing trend in downslope ndvi and no change in upslope ndvi  
downslope_decrease_upslope_nonsig <- downslope_trends[p_value <= 0.05 & slope < 0, ]$wsid %in% upslope_trends[p_value > 0.05,]$wsid
middle_dt <- downslope_trends[p_value <= 0.05 & slope < 0, ]
sum( middle_dt[downslope_decrease_upslope_nonsig,]$wsid %in% ratio_trends[p_value <= 0.05 & slope < 0, ]$wsid)/length(ratio_trends[p_value <= 0.05 & slope < 0, ]$wsid)  # 0.085
sum( middle_dt[downslope_decrease_upslope_nonsig,]$wsid %in% sd_trends[p_value <= 0.05 & slope < 0, ]$wsid)/length(sd_trends[p_value <= 0.05 & slope < 0, ]$wsid)  # 0.019

# make a csv for arcpro 
sig_trends <- ndvi_trends[p_value <= 0.05 & metric == "trend_ratio_ndvi",]
sig_trends[slope < 0,]$slope <- -1
sig_trends[slope > 0,]$slope <- 1
sig_trends <- sig_trends[,c("wsid", "slope")]
sig_trends$wsid <- as.numeric(sig_trends$wsid)
catchments <- st_read(paste0(home, "Data/Catchments/Headwater/headwater_catchments_perm_forest_32617.shp"))
catchments <- merge(catchments, sig_trends, by.x="NHDPlusID", by.y="wsid", inner=T)
st_write(catchments, paste0(home, "Data/Catchments/Headwater/headwater_catchments_sig_trend.shp"))

#################################################################################################################################################
#################################################################################################################################################
#################################################################################################################################################
# Create climate df
# Climate variables to include: tmin, tmax, vpd, par, prcp 

get_climate <- function(VAR){
  dt <- fread(paste0(home, "/Data/Climate/Summary/gs_", VAR, ".csv"))
  anomalies <- dt[!Type == "climatology", ]
  anomalies$Var_Name <- paste0(anomalies$Var, "_", anomalies$Season)
  anomalies <- anomalies[,c("NHDPlusID", "Value", "Year", "Var_Name")]
  anomalies_wide <- reshape(data=anomalies, idvar=c("NHDPlusID", "Year"), timevar="Var_Name", v.names="Value", direction="wide")
  new_column_names <- c(colnames(anomalies_wide)[1:2],substr(colnames(anomalies_wide)[3:length(colnames(anomalies_wide))], 7, nchar(colnames(anomalies_wide)[3:length(colnames(anomalies_wide))])))
  colnames(anomalies_wide) <- new_column_names
  return(anomalies_wide)
}

tmin <- get_climate("tmin")
tmax <- get_climate("tmax")
prcp <- get_climate("prcp")
vp <- get_climate("vp")
vpd <- get_climate("vpd")
par <- get_climate("par")
clim_dt <- Reduce(function(x, y) merge(x, y, all=TRUE, by=c("NHDPlusID", "Year")), list(tmin, tmax, prcp, vp, vpd, par))
rm(tmin, tmax, prcp, vp, vpd, par)

# get the R2 for each catchment of tmin and vp 
get_r2 <- function(X, Y){
  mod <- lm(Y~X)
  r2 <- summary(mod)$r.squared
  return(r2)
}

tmin_vp_r2 <- clim_dt[, .(R2=get_r2(tmin_AMJJA, vp_AMJJA)), .(NHDPlusID)]
mean(tmin_vp_r2$R2) # super high, so we drop vp 

#################################################################################################################################################
#################################################################################################################################################
#################################################################################################################################################
## Get the annual ndvi ratio, ndvi sd, ndvi downslope, ndvi upslope 
## NHDPlusID, Year, ratio
## had to do this part on hpc 
home <- "/share/klmarti3/kmcquil/Chapter2_mechanisms_forest_water_cycling"
ndvi_files <- list.files(paste0(home, "/Data/NDVI/catchment_ratio_ndvi_results"), full.names=TRUE, pattern="*.csv$")
ndvi_files <- unique(ndvi_files)

for(i in 1:length(ndvi_files)){
  if(file.size(ndvi_files[i]) == 0L){
    next
  }
  file <- fread(ndvi_files[i])
  file$year <- as.numeric(substr(file$Date, 1, 4))
  
  # 3 sigma rule 
  ratio_low <- mean(file$RatioDownUp_NDVI) - (3*sd(file$RatioDownUp_NDVI))
  ratio_high <- mean(file$RatioDownUp_NDVI) + (3*sd(file$RatioDownUp_NDVI))
  sd_low <- mean(file$SD_NDVI) - (3*sd(file$SD_NDVI))
  sd_high <- mean(file$SD_NDVI) + (3*sd(file$SD_NDVI))
  up_low <- mean(file$Upslope_NDVI) - (3*sd(file$Upslope_NDVI))
  up_high <- mean(file$Upslope_NDVI) + (3*sd(file$Upslope_NDVI))
  down_low <- mean(file$Downslope_NDVI) - (3*sd(file$Downslope_NDVI))
  down_high <- mean(file$Downslope_NDVI) + (3*sd(file$Downslope_NDVI))
  file <- file[!(RatioDownUp_NDVI < ratio_low | RatioDownUp_NDVI > ratio_high | SD_NDVI < sd_low | SD_NDVI > sd_high |
                   Upslope_NDVI < up_low | Upslope_NDVI > up_high | Downslope_NDVI < down_low | Downslope_NDVI > down_high),]
  
  summary <- file[,.(RatioDownUp_NDVI=mean(RatioDownUp_NDVI), SD_NDVI=mean(SD_NDVI), Downslope_NDVI=mean(Downslope_NDVI), Upslope_NDVI=mean(Upslope_NDVI),
                     Mean_NDVI = mean(Mean_NDVI), WSID=mean(WSID)), .(year)]
  if(i == 1){
    result <- summary
  }else{
    result <- rbind(result, summary)
  }
}

fwrite(result, paste0(home, "/Data/NDVI/annual_ratio_ndvi.csv"))

#################################################################################################################################################
#################################################################################################################################################
#################################################################################################################################################
# Create the dataframe of response and standardized predictors 
ndvi_dt <- fread(paste0(home, "/Data/NDVI/annual_ratio_ndvi.csv"))
colnames(ndvi_dt) <- c("Year", "RatioDownUp_NDVI", "SD_NDVI", "Downslope_NDVI", "Upslope_NDVI","Mean_NDVI","NHDPlusID")
ndvi_dt$NHDPlusID <- as.numeric(ndvi_dt$NHDPlusID)

ratio_ndvi_mean <- ndvi_dt[, .(mean_ratio_ndvi=mean(RatioDownUp_NDVI), 
                               sd_ratio_ndvi = sd(RatioDownUp_NDVI), 
                               mean_mean_ndvi = mean(Mean_NDVI), 
                               sd_mean_ndvi = sd(Mean_NDVI), 
                               mean_sd_ndvi=mean(SD_NDVI),
                               sd_sd_ndvi=sd(SD_NDVI)), .(NHDPlusID)]
ndvi_dt <- merge(ndvi_dt, ratio_ndvi_mean, by='NHDPlusID', left=TRUE)
ndvi_dt$RatioDownUp_NDVI_zscore <- (ndvi_dt$RatioDownUp_NDVI-ndvi_dt$mean_ratio_ndvi)/ndvi_dt$sd_ratio_ndvi
ndvi_dt$Mean_NDVI_zscore <- (ndvi_dt$Mean_NDVI - ndvi_dt$mean_mean_ndvi)/ndvi_dt$sd_mean_ndvi
ndvi_dt$SD_NDVI_zscore <- (ndvi_dt$SD_NDVI - ndvi_dt$mean_sd_ndvi)/ndvi_dt$sd_sd_ndvi

# merge with the climate variables
ndvi_cols_keep <- c("NHDPlusID", "Year", "RatioDownUp_NDVI", "Mean_NDVI", "RatioDownUp_NDVI_zscore", "SD_NDVI_zscore","Mean_NDVI_zscore")
full_dt <- merge(clim_dt, ndvi_dt[,..ndvi_cols_keep], by = c("NHDPlusID", "Year"), inner = TRUE)

# get the centroid of each polygon and join by NHDPlusID 
roi <- st_read(paste0(home, "/Data/Catchments/Headwater/headwater_catchments_perm_forest_32617.shp"))
roi <- cbind(roi, st_coordinates(st_centroid(roi)))
roi <- st_drop_geometry(roi[,c("NHDPlusID", "X", "Y")])
full_dt <- merge(full_dt, roi, by = "NHDPlusID")

#################################################################################################################################################
#################################################################################################################################################
#################################################################################################################################################
# Calculate the slope, r2, and p-value of each catchment NDVI with each climate variable and save

vars <- c("tmin", "tmax", "prcp", "vpd", "par")
seasons <- c("AMJJA")
ids <- unique(full_dt$NHDPlusID)

K <- 0
for(var in vars){
  for(season in seasons){
    clim_column<- paste0(var, "_", season)
    for(id in ids){
      K <- K+1
      corr <- lm(full_dt[NHDPlusID == id,get("RatioDownUp_NDVI_zscore")] ~ full_dt[NHDPlusID == id,get(clim_column)])
      if(K == 1){
        results <- data.table('NHDPlusID'=id, 
                              'Season'=season,
                              'Variable'=var,
                              'Slope'=corr$coefficients[2],
                              'Pvalue'=summary(corr)$coefficients[2,4], 
                              'Rsquared'=summary(corr)$r.squared)
      }else{
        row = data.table('NHDPlusID'=id, 
                         'Season'=season,
                         'Variable'=var,
                         'Slope'=corr$coefficients[2],
                         'Pvalue'=summary(corr)$coefficients[2,4], 
                         'Rsquared'=summary(corr)$r.squared)
        results <- rbind(results, row)
      }
    }
  }
}
fwrite(results, paste0(home, "/Data/growing_season_climate_NDVI_linmod_r2.csv"))

K <- 0
for(var in vars){
  for(season in seasons){
    clim_column<- paste0(var, "_", season)
    for(id in ids){
      K <- K+1
      corr <- lm(full_dt[NHDPlusID == id,get("SD_NDVI_zscore")] ~ full_dt[NHDPlusID == id,get(clim_column)])
      if(K == 1){
        results <- data.table('NHDPlusID'=id, 
                              'Season'=season,
                              'Variable'=var,
                              'Slope'=corr$coefficients[2],
                              'Pvalue'=summary(corr)$coefficients[2,4], 
                              'Rsquared'=summary(corr)$r.squared)
      }else{
        row = data.table('NHDPlusID'=id, 
                         'Season'=season,
                         'Variable'=var,
                         'Slope'=corr$coefficients[2],
                         'Pvalue'=summary(corr)$coefficients[2,4], 
                         'Rsquared'=summary(corr)$r.squared)
        results <- rbind(results, row)
      }
    }
  }
}
fwrite(results, paste0(home, "/Data/growing_season_climate_SD_NDVI_linmod_r2.csv"))


results <- fread(paste0(home, "/Data/growing_season_climate_NDVI_linmod_r2.csv"))
#results <- fread(paste0(home, "/Data/growing_season_climate_SD_NDVI_linmod_r2.csv"))
results$Group <- paste0(results$Variable, " ", results$Season)
results$Group <- as.factor(results$Group)
results_sig <- results[results$Pvalue <= 0.05,]
summary_results_sig <- results_sig[, .(count = .N/30044), .(Group)]

# with the gridmet data vpd was significant in 7% of catchments 
ggplot(results_sig, aes(x=Rsquared)) + 
  geom_histogram(binwidth=0.1, color='sky blue', fill='sky blue') + 
  geom_vline(xintercept = 0, color='black', linetype='dashed') + 
  theme_bw() + 
  facet_wrap(~Group)

ggplot(results_sig, aes(x=Slope)) + 
  geom_histogram(binwidth=0.1, color='sky blue', fill='sky blue') + 
  geom_vline(xintercept = 0, color='black', linetype='dashed') + 
  theme_bw() + 
  facet_wrap(~Group)


#################################################################################################################################################
#################################################################################################################################################
#################################################################################################################################################
# create a figure mapping the overall ratio NDVI and the trend in ratio NDVI 
catch <- st_read(paste0(home, "/Data/Catchments/Headwater/headwater_catchments_perm_forest_32617.shp"))
catch <- catch[,c("NHDPlusID", "geometry")]
roi <- st_read("/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling/Data/ROI/blue_ridge_plus_reference.shp")
roi <- st_transform(roi, crs="EPSG:32617")

ndvi_dt <- fread(paste0(home, "/Data/NDVI/annual_ratio_ndvi.csv"))
ndvi_summary <- ndvi_dt[,.(RatioDownUp_NDVI = mean(RatioDownUp_NDVI)), .(WSID)]
ndvi_summary$WSID <- as.numeric(ndvi_summary$WSID)
catch_summary <- merge(catch, ndvi_summary, by.x="NHDPlusID", by.y="WSID", inner=TRUE)

nrow(catch_summary[catch_summary$RatioDownUp_NDVI < 1,])/nrow(catch_summary)

grp_cors <- cbind(catch_summary, st_coordinates(st_centroid(catch_summary)))
grp_cors$Group <- rep(" ", nrow(grp_cors))
ratio_violin_plot <- ggplot(grp_cors, aes(x=Group,y=RatioDownUp_NDVI)) +
  geom_violin() + 
  coord_flip() + 
  theme_minimal() + 
  xlab("") + 
  ylab("") + 
  theme(axis.title.y = element_text(size = rel(0.1)))

upper_limit <- quantile(grp_cors$RatioDownUp_NDVI, 0.975)
lower_limit <- quantile(grp_cors$RatioDownUp_NDVI, 0.025)
grp_cors[grp_cors$RatioDownUp_NDVI < lower_limit,]$RatioDownUp_NDVI <- lower_limit
grp_cors[grp_cors$RatioDownUp_NDVI > upper_limit,]$RatioDownUp_NDVI <- upper_limit
ratio_map <- ggplot() + 
  geom_point(data=grp_cors, aes(x=X, y=Y, color=RatioDownUp_NDVI), alpha=0.5, size = 0.05) + 
  #geom_sf(data=grp_cors, aes(fill=RatioDownUp_NDVI, color=RatioDownUp_NDVI),size=0.1) + 
  scale_color_gradient2(low = "#F62C00", mid = "light yellow", high = "#440154", midpoint = 1, name="", limits=c(lower_limit, upper_limit)) + 
  scale_fill_gradient2(low = "#F62C00", mid = "light yellow", high = "#440154", midpoint = 1, name="", limits=c(lower_limit, upper_limit)) + 
  geom_sf(data=roi, fill=NA, color='black')+ 
  theme_classic() + 
  ggtitle(bquote(NDVI[Down:Up])) + 
  xlab("") + 
  ylab("") + 
  theme(legend.position = c(0.9, 0.28)) + 
  inset_element(ratio_violin_plot, 0.4, 0.07, 0.75, 0.35, align_to = 'full') #left bottom right top

# now make the trend ratio map 
trend_dt <- fread(paste0(home, "/Data/NDVI/ratio_ndvi_trend_results_outliers.csv"))
#lower_limit <- quantile(trend_dt[metric %in% c("trend_ratio_ndvi", "trend_upslope_ndvi", "trend_downslope_ndvi"),]$slope, 0.025)
#upper_limit <- quantile(trend_dt[metric %in% c("trend_ratio_ndvi", "trend_upslope_ndvi", "trend_downslope_ndvi"),]$slope, 0.975)
lower_limit <- quantile(trend_dt[metric %in% c("trend_ratio_ndvi"),]$slope, 0.025)
upper_limit <- quantile(trend_dt[metric %in% c("trend_ratio_ndvi"),]$slope, 0.975)

resp <- trend_dt[trend_dt$metric == "trend_ratio_ndvi",]
resp <- resp[resp$p_value <= 0.05,]
resp <- resp[,c("wsid", "slope")]
resp$wsid <- as.numeric(resp$wsid)
resp_summary <- merge(catch, resp, by.x="NHDPlusID", by.y="wsid", inner=TRUE)
resp_sf <- cbind(resp_summary, st_coordinates(st_centroid(resp_summary)))

resp_sf$Group <- rep(" ", nrow(resp_sf))
trend_violin_plot <- ggplot(resp_sf, aes(x=Group,y=slope)) +
  geom_violin() + 
  coord_flip() + 
  theme_minimal() + 
  xlab("") + 
  ylab("") + 
  theme(axis.title.y = element_text(size = rel(0.1)))

resp_sf[resp_sf$slope < lower_limit,]$slope <- lower_limit
resp_sf[resp_sf$slope > upper_limit,]$slope <- upper_limit
trend_map <- ggplot() + 
  geom_point(data=resp_sf, aes(x=X, y=Y, color=slope), alpha=0.7, size = 0.05) + 
  #geom_sf(data=resp_sf, aes(fill=slope, color=slope), size=0.1) + 
  #geom_point(data=grp_sig, aes(x=X, y=Y, color=Correlation, size = abs(Correlation)), alpha=0.4) + 
  scale_color_gradient2(low = "#F62C00", mid = "light yellow", high = "#440154", midpoint = 0, name="", limits=c(lower_limit, upper_limit)) + 
  scale_fill_gradient2(low = "#F62C00", mid = "light yellow", high = "#440154", midpoint = 0, name="", limits=c(lower_limit, upper_limit)) + 
  geom_sf(data=roi, fill=NA, color='black')+ 
  #scale_size(range=c(0, 5)) + 
  theme_classic() + 
  ggtitle(bquote('Trend ' * NDVI[Down:Up])) + 
  #ggtitle("") + 
  xlab("") + 
  ylab("") + 
  #theme(legend.position = "none") + 
  theme(legend.position = c(0.9, 0.28)) + 
  inset_element(trend_violin_plot, 0.4, 0.07, 0.8, 0.35, align_to = 'full') #left bottom right top


tiff(paste0(home, "/Figures/ratio_trend_ndvi_maps.tiff"), units = 'in', res=1000, width = 9, height = 3.75, compression="lzw")
ggarrange(ratio_map, trend_map, nrow = 1, ncol = 2, common.legend = FALSE, labels="AUTO")
dev.off()

tiff(paste0(home, "/Figures/Final_Figures/Fig3.tiff"), units = 'in', res=1000, width = 9, height = 3.75, compression="lzw")
ggarrange(ratio_map, trend_map, nrow = 1, ncol = 2, common.legend = FALSE, labels="AUTO")
dev.off()
# manually go into preview on mac and adjust size 

# now make the trend SD map 
trend_dt <- fread(paste0(home, "/Data/NDVI/ratio_ndvi_trend_results_outliers.csv"))
lower_limit <- quantile(trend_dt[metric %in% c("trend_sd_ndvi"),]$slope, 0.025)
upper_limit <- quantile(trend_dt[metric %in% c("trend_sd_ndvi"),]$slope, 0.975)

resp <- trend_dt[trend_dt$metric == "trend_sd_ndvi",]
resp <- resp[resp$p_value <= 0.05,]
resp <- resp[,c("wsid", "slope")]
resp$wsid <- as.numeric(resp$wsid)
resp_summary <- merge(catch, resp, by.x="NHDPlusID", by.y="wsid", inner=TRUE)
resp_sf <- cbind(resp_summary, st_coordinates(st_centroid(resp_summary)))

resp_sf$Group <- rep(" ", nrow(resp_sf))
trend_violin_plot <- ggplot(resp_sf, aes(x=Group,y=slope)) +
  geom_violin() + 
  coord_flip() + 
  theme_minimal() + 
  xlab("") + 
  ylab("") + 
  theme(axis.title.y = element_text(size = rel(0.1))) + 
  theme(text = element_text(size = 7))     

resp_sf[resp_sf$slope < lower_limit,]$slope <- lower_limit
resp_sf[resp_sf$slope > upper_limit,]$slope <- upper_limit
trend_map <- ggplot() + 
  geom_sf(data=roi, fill='white', color='black')+ 
  geom_point(data=resp_sf, aes(x=X, y=Y, color=slope), alpha=0.7, size = 0.05) + 
  scale_color_gradient2(low = "#F62C00", mid = "light yellow", high = "#440154", midpoint = 0, name="", limits=c(lower_limit, upper_limit)) + 
  theme_classic() + 
  ggtitle('Trend SD NDVI') + 
  xlab("") + 
  ylab("") + 
  theme(legend.position = c(0.9, 0.28)) + 
  inset_element(trend_violin_plot, 0.4, 0.07, 0.80, 0.35, align_to = 'full') #left bottom right top

tiff(paste0(home, "/Figures/SD_trend_ndvi_maps.tiff"), units = 'in', res=300, width = 5, height = 3.75)
trend_map
dev.off()



# make the downslope trend map 
lower_limit <- quantile(trend_dt[metric %in% c("trend_upslope_ndvi", "trend_downslope_ndvi"),]$slope, 0.025)
upper_limit <- quantile(trend_dt[metric %in% c("trend_upslope_ndvi", "trend_downslope_ndvi"),]$slope, 0.975)

downslope <- trend_dt[metric == "trend_downslope_ndvi",]
downslope <- downslope[p_value <= 0.05,]
downslope <- downslope[,c("wsid", "slope")]
downslope$wsid <- as.numeric(downslope$wsid)
downslope_summary <- merge(catch, downslope, by.x="NHDPlusID", by.y='wsid', inner=TRUE)
downslope_sf <- cbind(downslope_summary, st_coordinates(st_centroid(downslope_summary)))
downslope_sf$Group <- rep(" ", nrow(downslope_sf))
downslope_violin_plot <- ggplot(downslope_sf, aes(x=Group,y=slope)) +
  geom_violin() + 
  coord_flip() + 
  theme_minimal() + 
  xlab("") + 
  ylab("") + 
  theme(axis.title.y = element_text(size = rel(0.1)))
#lower_limit <- quantile(downslope_sf$slope, 0.05)
#upper_limit <- quantile(downslope_sf$slope, 0.95)
downslope_sf[downslope_sf$slope < lower_limit,]$slope <- lower_limit
downslope_sf[downslope_sf$slope > upper_limit,]$slope <- upper_limit
downslope_map <- ggplot() + 
  geom_sf(data=roi, fill='white', color='black')+ 
  geom_point(data=downslope_sf, aes(x=X, y=Y, color=slope), alpha=0.7, size = 0.05) + 
  #geom_point(data=grp_sig, aes(x=X, y=Y, color=Correlation, size = abs(Correlation)), alpha=0.4) + 
  scale_color_gradient2(low = "#F62C00", mid = "light yellow", high = "#440154", midpoint = 0, name="", limits=c(lower_limit, upper_limit)) + 
  #scale_size(range=c(0, 5)) + 
  theme_classic() + 
  ggtitle('Trend Downslope NDVI') + 
  xlab("") + 
  ylab("") + 
  #theme(legend.position = c(0.85, 0.3)) + 
  theme(legend.position = "none") + 
  inset_element(downslope_violin_plot, 0.4, 0.07, 0.95, 0.35, align_to = 'full') #left bottom right top

# make the downslope trend map 
upslope <- trend_dt[metric == "trend_upslope_ndvi",]
upslope <- upslope[p_value <= 0.05,]
upslope <- upslope[,c("wsid", "slope")]
upslope$wsid <- as.numeric(upslope$wsid)
upslope_summary <- merge(catch, upslope, by.x="NHDPlusID", by.y='wsid', inner=TRUE)
upslope_sf <- cbind(upslope_summary, st_coordinates(st_centroid(upslope_summary)))
upslope_sf$Group <- rep(" ", nrow(upslope_sf))
upslope_violin_plot <- ggplot(upslope_sf, aes(x=Group,y=slope)) +
  geom_violin() + 
  coord_flip() + 
  theme_minimal() + 
  xlab("") + 
  ylab("") + 
  theme(axis.title.y = element_text(size = rel(0.1)))
#lower_limit <- quantile(upslope_sf$slope, 0.05)
#upper_limit <- quantile(upslope_sf$slope, 0.95)
upslope_sf[upslope_sf$slope < lower_limit,]$slope <- lower_limit
upslope_sf[upslope_sf$slope > upper_limit,]$slope <- upper_limit
upslope_map <- ggplot() + 
  geom_sf(data=roi, fill='white', color='black')+ 
  geom_point(data=upslope_sf, aes(x=X, y=Y, color=slope), alpha=0.7, size = 0.05) + 
  #geom_point(data=grp_sig, aes(x=X, y=Y, color=Correlation, size = abs(Correlation)), alpha=0.4) + 
  scale_color_gradient2(low = "#F62C00", mid = "light yellow", high = "#440154", midpoint = 0, name='', limits=c(lower_limit, upper_limit)) + 
  #scale_size(range=c(0, 5)) + 
  theme_classic() + 
  ggtitle('Trend Upslope NDVI') + 
  xlab("") + 
  ylab("") + 
  theme(legend.position = c(0.9, 0.3)) +
  #theme(legend.position ="none") +
  inset_element(upslope_violin_plot, 0.4, 0.07, 0.75, 0.35, align_to = 'full') #left bottom right top

#tiff(paste0(home, "/Figures/ratio_trend_ndvi_maps.tiff"), units = 'in', res=300, width = 9, height = 7.5)
#tiff(paste0(home, "/Figures/updown_trend_ndvi_maps.tiff"), units = 'in', res=300, width = 9, height = 3.75)
tiff(paste0(home, "/Figures/updown_trend_ndvi_maps.tiff"), units = 'in', res=300, width = 9, height = 7.5)
#ggarrange(ratio_map, trend_map,downslope_map, upslope_map, nrow = 2, ncol = 2, common.legend = TRUE, legend.grob=legend, legend='right')
#ggarrange(ratio_map, trend_map,downslope_map, upslope_map, nrow = 2, ncol = 2, common.legend = FALSE, labels="AUTO")
#ggarrange(downslope_map, upslope_map, nrow = 2, ncol = 2, common.legend = FALSE, labels="AUTO")
ggarrange(downslope_map, upslope_map, trend_map,nrow = 2, ncol = 2, common.legend = FALSE, labels="AUTO")
dev.off()


#################################################################################################################################################
#################################################################################################################################################
#################################################################################################################################################
# Create a bar chart of the fraction of catchments with a decreasing, increasing, non sig trend for ratio, SD, downslope, uslope 

ratio_frac_dt <- data.table(var=c("Down:Up", "Down:Up", "Down:Up"), 
                            group = c("Increasing", "Decreasing", "No Change"), 
                            value = c(nrow(trend_dt[trend_dt$metric == "trend_ratio_ndvi" & slope > 0 & p_value<=0.05,])/nrow(trend_dt[trend_dt$metric == "trend_ratio_ndvi",]), 
                                      nrow(trend_dt[trend_dt$metric == "trend_ratio_ndvi" & slope < 0 & p_value<=0.05,])/nrow(trend_dt[trend_dt$metric == "trend_ratio_ndvi",]), 
                                      nrow(trend_dt[trend_dt$metric == "trend_ratio_ndvi" & p_value>0.05,])/nrow(trend_dt[trend_dt$metric == "trend_ratio_ndvi",])))

sd_frac_dt <- data.table(var=c("SD", "SD", "SD"), 
                            group = c("Increasing", "Decreasing", "No Change"), 
                            value = c(nrow(trend_dt[trend_dt$metric == "trend_sd_ndvi" & slope > 0 & p_value<=0.05,])/nrow(trend_dt[trend_dt$metric == "trend_sd_ndvi",]), 
                                      nrow(trend_dt[trend_dt$metric == "trend_sd_ndvi" & slope < 0 & p_value<=0.05,])/nrow(trend_dt[trend_dt$metric == "trend_sd_ndvi",]), 
                                      nrow(trend_dt[trend_dt$metric == "trend_sd_ndvi" & p_value>0.05,])/nrow(trend_dt[trend_dt$metric == "trend_sd_ndvi",])))

upslope_frac_dt <- data.table(var=c("Upslope", "Upslope", "Upslope"), 
                         group = c("Increasing", "Decreasing", "No Change"), 
                         value = c(nrow(trend_dt[trend_dt$metric == "trend_upslope_ndvi" & slope > 0 & p_value<=0.05,])/nrow(trend_dt[trend_dt$metric == "trend_upslope_ndvi",]), 
                                   nrow(trend_dt[trend_dt$metric == "trend_upslope_ndvi" & slope < 0 & p_value<=0.05,])/nrow(trend_dt[trend_dt$metric == "trend_upslope_ndvi",]), 
                                   nrow(trend_dt[trend_dt$metric == "trend_upslope_ndvi" & p_value>0.05,])/nrow(trend_dt[trend_dt$metric == "trend_upslope_ndvi",])))

downslope_frac_dt <- data.table(var=c("Downslope", "Downslope", "Downslope"), 
                              group = c("Increasing", "Decreasing", "No Change"), 
                              value = c(nrow(trend_dt[trend_dt$metric == "trend_downslope_ndvi" & slope > 0 & p_value<=0.05,])/nrow(trend_dt[trend_dt$metric == "trend_downslope_ndvi",]), 
                                        nrow(trend_dt[trend_dt$metric == "trend_downslope_ndvi" & slope < 0 & p_value<=0.05,])/nrow(trend_dt[trend_dt$metric == "trend_downslope_ndvi",]), 
                                        nrow(trend_dt[trend_dt$metric == "trend_downslope_ndvi" & p_value>0.05,])/nrow(trend_dt[trend_dt$metric == "trend_downslope_ndvi",])))

frac_dt <- rbind(ratio_frac_dt, sd_frac_dt, upslope_frac_dt, downslope_frac_dt)
frac_dt$var <- factor(frac_dt$var, levels=c("Down:Up", "SD", "Downslope", "Upslope"))
frac_dt$value <- frac_dt$value * 100

frac_plot <- ggplot(frac_dt, aes(x=var, y=value, fill=group)) + 
  geom_bar(stat='identity', position='dodge') + 
  scale_fill_manual(values=c("#7570b3", "#bf5b17", "#a6cee3"), name="Trend") + 
  theme_classic() + 
  xlab("NDVI") + 
  ylab("% Catchments") + 
  theme(legend.position = c(0.25, 0.85))

tiff(paste0(home, "/Figures/fraction_ndvi_trends.tiff"), units = 'in', res=1000, width = 5, height = 4)
frac_plot
dev.off()

tiff(paste0(home, "/Figures/Final_Figures/Fig4.tiff"), units = 'in', res=1000, width = 5, height = 4)
frac_plot
dev.off()
#manually resize in preview on mac 

#################################################################################################################################################
#################################################################################################################################################
#################################################################################################################################################
# create a figure mapping the slope and rsquared for each catchment x gs climate variable
# include a violin plot with each map to get a better sense of the distribution 
# only include catchments with significant relationships 

# start with the Slope 
# get the upper and lower bounds for the color bar 
upper_limit <- quantile(results_sig$Slope, 0.975)
lower_limit <- quantile(results_sig$Slope, 0.025)
plot_ndvi_clim_slopes <- function(grp, var, title){
  grp_cors <- merge(catch, results_sig[Group == grp], inner=TRUE, by="NHDPlusID")
  grp_cors <- cbind(grp_cors, st_coordinates(st_centroid(grp_cors)))
  
  grp_cors_vi <- grp_cors
  grp_cors_vi$Group <- rep("", nrow(grp_cors_vi))
  trend_violin_plot <- ggplot(grp_cors_vi, aes(x=Group,y=get(var))) +
    geom_violin() + 
    coord_flip() + 
    theme_minimal() + 
    xlab("") + 
    ylab("") + 
    theme(axis.title.y = element_text(size = rel(0.1)))
  
  grp_cors[grp_cors$Slope < lower_limit, "Slope"] <- lower_limit
  grp_cors[grp_cors$Slope > upper_limit, "Slope"] <- upper_limit
  
  gp <- ggplot() + 
    geom_sf(data=roi, fill='white', color='black')+ 
    geom_point(data=grp_cors, aes(x=X, y=Y, color=get(var)), alpha=0.4, size=0.01) + 
    theme_classic() + 
    xlab("") + 
    ylab("") + 
    ggtitle(title) +
    #scale_color_gradient2(low = "#F62C00", mid = "light yellow", high = "#440154", midpoint = 0, name=var, limits=c(lower_limit, upper_limit)) +
    #theme(legend.title = element_blank()) + 
    #theme(legend.position = "none") +
    #inset_element(trend_violin_plot, 0.5, 0.07, 1, 0.40, align_to = 'full') #left bottom right top
    scale_color_gradient2(low = "#F62C00", mid = "light yellow", high = "#440154", midpoint = 0, name="Slope", limits=c(lower_limit, upper_limit)) +
    theme(legend.position = c(0.95, 0.3)) + 
    inset_element(trend_violin_plot, 0.45, 0.07, 0.9, 0.40, align_to = 'full') #left bottom right top
  return(gp)
}

results_sig$NHDPlusID <- as.numeric(results_sig$NHDPlusID)
catch <- st_read(paste0(home, "/Data/Catchments/Headwater/headwater_catchments_perm_forest_32617.shp"))
catch <- catch[,c("NHDPlusID", "geometry")]
roi <- st_read("/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling/Data/ROI/blue_ridge_plus_reference.shp")
roi <- st_transform(roi, crs="EPSG:32617")

tmin_slope <- plot_ndvi_clim_slopes("tmin AMJJA", "Slope", bquote(T[min]))
tmax_slope <- plot_ndvi_clim_slopes("tmax AMJJA", "Slope", bquote(T[max]))
vpd_slope <- plot_ndvi_clim_slopes("vpd AMJJA", "Slope", "VPD")
prcp_slope <- plot_ndvi_clim_slopes("prcp AMJJA", "Slope", "Prcp")
par_slope <- plot_ndvi_clim_slopes("par AMJJA", "Slope", "PAR")

amjja_results_sig <- results_sig[grep("par AMJJA", results_sig$Group),]
amjja_results_sig <- merge(catch, amjja_results_sig, inner=TRUE, by="NHDPlusID")
amjja_results_sig <- cbind(amjja_results_sig, st_coordinates(st_centroid(amjja_results_sig)))
legend <- ggplot() + 
  geom_sf(data=roi, fill='white', color='black')+ 
  geom_point(data=amjja_results_sig, aes(x=X, y=Y, color=Slope), alpha=0.4, size=0.01) + 
  scale_color_gradient2(low = "#F62C00", mid = "light yellow", high = "#440154", midpoint = 0, name='Slope', limits=c(lower_limit, upper_limit)) +
  theme_classic() + 
  xlab("") + 
  ylab("") + 
  ggtitle("") + 
  #theme(legend.title = element_blank()) + 
  theme(legend.position = "right", 
        legend.text=element_text(size=12)) +
  theme(axis.text=element_text(size=12, color = 'black'),
        axis.title=element_text(size=12,face="bold")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1)) 
legend <- ggpubr::get_legend(legend)

tiff(paste0(home, "/Figures/slope_ndvi_clim_models.tiff"), units = 'in', res=1000, width = 10, height = 12)
#ggarrange(tmin_slope, tmax_slope, vpd_slope, prcp_slope, par_slope, legend, nrow = 3, ncol = 2, labels='AUTO') 
ggarrange(tmin_slope, tmax_slope, vpd_slope, prcp_slope, par_slope, nrow = 3, ncol = 2, labels='AUTO') 
          #common.legend = TRUE,
          #legend.grob=legend,
          #legend = "right", 
          #labels='AUTO')
#ggarrange(tmin_slope, tmax_slope, vpd_slope, prcp_slope, par_slope, nrow = 3, ncol = 2, 
#          common.legend = TRUE,
#          legend.grob=legend,
#          legend = "right", 
#          labels='AUTO')
dev.off()

tiff(paste0(home, "/Figures/Final_Figures/Fig6.tiff"), units = 'in', res=1000, width = 10, height = 12)
ggarrange(tmin_slope, tmax_slope, vpd_slope, prcp_slope, par_slope, nrow = 3, ncol = 2, labels='AUTO') 
dev.off()
# resize on mac in preview 


# now create the figure for the R-squared 
upper_limit <- quantile(results_sig$Rsquared, 0.975)
lower_limit <- quantile(results_sig$Rsquared, 0.025)
plot_ndvi_clim_rsquared <- function(grp, var, title){
  grp_cors <- merge(catch, results_sig[Group == grp], inner=TRUE, by="NHDPlusID")
  grp_cors <- cbind(grp_cors, st_coordinates(st_centroid(grp_cors)))
  
  grp_cors_vi <- grp_cors
  grp_cors_vi$Group <- rep("", nrow(grp_cors_vi))
  trend_violin_plot <- ggplot(grp_cors_vi, aes(x=Group,y=get(var))) +
    geom_violin() + 
    coord_flip() + 
    theme_minimal() + 
    xlab("") + 
    ylab("") + 
    theme(axis.title.y = element_text(size = rel(0.1)))
  
  grp_cors[grp_cors$Rsquared < lower_limit, "Rsquared"] <- lower_limit
  grp_cors[grp_cors$Rsquared > upper_limit, "Rsquared"] <- upper_limit
  
  gp <- ggplot() + 
    geom_sf(data=roi, fill='white', color='black')+ 
    geom_point(data=grp_cors, aes(x=X, y=Y, color=get(var)), alpha=0.4, size=0.01) + 
    theme_classic() + 
    xlab("") + 
    ylab("") + 
    ggtitle(title) +
    theme(legend.title = element_blank()) + 
    theme(legend.position = "none") +
    scale_color_viridis_c(limits = c(lower_limit, upper_limit)) + 
    inset_element(trend_violin_plot, 0.5, 0.07, 1, 0.40, align_to = 'full')   #left bottom right top
    #scale_color_viridis_c(limits = c(lower_limit, upper_limit),name=bquote(R^2)) + 
    #theme(legend.position = c(0.95, 0.3)) + 
    #inset_element(trend_violin_plot, 0.45, 0.07, 0.9, 0.40, align_to = 'full') #left bottom right top
  return(gp)
}


tmin_rs <- plot_ndvi_clim_rsquared("tmin AMJJA", "Rsquared", bquote(T[min]))
tmax_rs <- plot_ndvi_clim_rsquared("tmax AMJJA", "Rsquared", bquote(T[max]))
vpd_rs <- plot_ndvi_clim_rsquared("vpd AMJJA", "Rsquared", "VPD")
prcp_rs <- plot_ndvi_clim_rsquared("prcp AMJJA", "Rsquared", "Precipitation")
par_rs <- plot_ndvi_clim_rsquared("par AMJJA", "Rsquared", "PAR")

legend <- ggplot() + 
  geom_sf(data=roi, fill='white', color='black')+ 
  geom_point(data=amjja_results_sig, aes(x=X, y=Y, color=Rsquared), alpha=0.4, size=0.01) + 
  scale_color_viridis_c(limits = c(lower_limit, upper_limit), name='Rsquared') + 
  theme_classic() + 
  xlab("") + 
  ylab("") + 
  ggtitle("") + 
  #theme(legend.title = element_blank()) + 
  theme(legend.position = "right", 
        legend.text=element_text(size=12)) +
  theme(axis.text=element_text(size=12, color = 'black'),
        axis.title=element_text(size=12,face="bold")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1)) 
legend <- ggpubr::get_legend(legend)

tiff(paste0(home, "/Figures/rsquared_ndvi_clim_models.tiff"), units = 'in', res=1000, width = 10, height = 12)
ggarrange(tmin_rs, tmax_rs, vpd_rs, prcp_rs, par_rs, nrow = 3, ncol = 2, labels='AUTO')
dev.off()

tiff(paste0(home, "/Figures/Final_Figures/Fig5.tiff"), units = 'in', res=1000, width = 10, height = 12)
ggarrange(tmin_rs, tmax_rs, vpd_rs, prcp_rs, par_rs, nrow = 3, ncol = 2, labels='AUTO')
dev.off()

#################################################################################################################################################
#################################################################################################################################################
#################################################################################################################################################
## Make a table of the quantiles of slopes, quantiles of rsquared, and the % of catchments with sig relationshihp 
## for each growing season climate variable 

#results <- fread(paste0(home, "/Data/growing_season_climate_NDVI_linmod_r2.csv"))
results <- fread(paste0(home, "/Data/growing_season_climate_SD_NDVI_linmod_r2.csv"))
results$Group <- paste0(results$Variable, " ", results$Season)
results$Group <- as.factor(results$Group)
results_sig <- results[results$Pvalue <= 0.05,]


amjja_results_sig_dt_summary <- results_sig[, .(fraction_sig = .N/30044, 
                                                   s1q=quantile(Slope, 0.25), 
                                                   smed=quantile(Slope, 0.5),
                                                   smean=mean(Slope),
                                                   s3q=quantile(Slope, 0.75), 
                                                   ssd=sd(Slope),
                                                   r1q=quantile(Rsquared, 0.25),
                                                   rmed=quantile(Rsquared, 0.5),
                                                   rmeaen=mean(Rsquared),
                                                   r3q=quantile(Rsquared, 0.75),
                                                   rsd=sd(Rsquared)), .(Group)]
amjja_results_sig_dt_summary[,2:ncol(amjja_results_sig_dt_summary)] <- round(amjja_results_sig_dt_summary[,2:ncol(amjja_results_sig_dt_summary)], 3)
fwrite(amjja_results_sig_dt_summary, paste0(home, "/Data/clim_SD_ndvi_model_summary.csv"))






#################################################################################################################################################
#################################################################################################################################################
#################################################################################################################################################
## Calculate the correlation between watershed characteristics and the trend in ratio downslope:upslope NDVI and overall downslope:upslope NDVI 

# climatology variables 
get_climatology <- function(VAR){
  dt <- fread(paste0(home, "/Data/Climate/Summary/gs_", VAR,".csv"))
  climatology <- dt[Type == "climatology" & Season == "AMJJA", c("NHDPlusID", "Value")]
  colnames(climatology) <- c("NHDPlusID", VAR)
  return(climatology)
}

tmin <- get_climatology("tmin")
tmax <- get_climatology("tmax")
prcp <- get_climatology("prcp")
vpd <- get_climatology("vpd")
par <- get_climatology("par")
climatology_dt <- Reduce(function(x, y) merge(x, y, all=TRUE), list(tmin, tmax, prcp, vpd, par))

# topography variables 
get_topo <- function(VAR){
  dt <- fread(paste0(home, "/Data/Topography/", VAR,"_hw_summary.csv"))
  colnames(dt) <- c("V1", "NHDPlusID", VAR)
  return(dt <- dt[,c(2,3)])
}
topo_dt <- Reduce(function(x,y) merge(x,y, all=TRUE), list(get_topo("Aspect"),get_topo("Slope"), get_topo("Elevation"), get_topo("Latitude"), get_topo("HLI")))

# % evergreen forest 
nlcd_2001 <- fread(paste0(home, "/Data/nlcd_permanent_forest/hw_nlcd_2001_summary.csv"))
nlcd_2016 <- fread(paste0(home, "/Data/nlcd_permanent_forest/hw_nlcd_2016_summary.csv"))

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

# get the area of each catchment
roi <- st_read(paste0(home, "/Data/Catchments/Headwater/headwater_catchments_perm_forest_32617.shp"))
area <- as.data.frame(roi)[,c('NHDPlusID', 'AreaSqKm_x')]

# put them all together
pred <- Reduce(function(x,y) merge(x,y,left=TRUE), list(climatology_dt, topo_dt, nlcd_summary[,c("NHDPlusID", "pct_evergreen")], area))

# Get the resposne variables (NDVI ratio, trend NDVI ratio, and mean NDVI)
resp <- fread(paste0(home, "/Data/NDVI/ratio_ndvi_trend_results_outliers.csv"))
resp <- resp[resp$metric == "trend_ratio_ndvi",]
resp$wsid <- as.numeric(resp$wsid)
resp_sub <- resp[,c("wsid", "slope", "mean_ndvi", "p_value")]
colnames(resp_sub) <- c("NHDPlusID", "Trend", "Mean_NDVI", "Pvalue_Trend")

colnames(ndvi_dt) <- c("year", "RatioDownUp_NDVI", "SD_NDVI", "Downslope_NDVI", "Upslope_NDVI","Mean_NDVI", "NHDPlusID")
overall_ndvi <- ndvi_dt[, .(Overall_RatioDownUp_NDVI=mean(RatioDownUp_NDVI)), .(NHDPlusID)]
overall_ndvi$NHDPlusID <- as.numeric(overall_ndvi$NHDPlusID)
resp_sub <- merge(resp_sub, overall_ndvi, by="NHDPlusID")

# combine the response and predictor variables and zscore standardize 
catchment_var_dt <- merge(resp_sub, pred)
catchment_var_dt <- cbind(catchment_var_dt[,c("NHDPlusID")], scale(catchment_var_dt[,2:ncol(catchment_var_dt)], center=TRUE, scale=TRUE))
fwrite(catchment_var_dt, paste0(home, "/Data/avg_catchment_vars.csv"))
catchment_var_dt <- catchment_var_dt[Pvalue_Trend <= 0.05,]


ugh <- function(x){
  hmm = ndvi_dt[ndvi_dt$year <= x, .(rdn = mean(RatioDownUp_NDVI)), .(NHDPlusID)]
  return(nrow(hmm[hmm$rdn >1, ])/nrow(hmm))
  
}
testing <- c()
for(i in 1986:2021){
  testing <- c(testing, ugh(i))
}


ratio_row <- data.table(
  Predictor = 'Downslope:Upslope NDVI Ratio', 
  Tmin = cor.test(catchment_var_dt$Overall_RatioDownUp_NDVI, catchment_var_dt$tmin, method = 'pearson')$estimate, 
  Tmax = cor.test(catchment_var_dt$Overall_RatioDownUp_NDVI, catchment_var_dt$tmax, method = 'pearson')$estimate, 
  VPD = cor.test(catchment_var_dt$Overall_RatioDownUp_NDVI, catchment_var_dt$vpd, method = 'pearson')$estimate, 
  Precipitation = cor.test(catchment_var_dt$Overall_RatioDownUp_NDVI, catchment_var_dt$prcp, method = 'pearson')$estimate, 
  PAR = cor.test(catchment_var_dt$Overall_RatioDownUp_NDVI, catchment_var_dt$par, method = 'pearson')$estimate,
  Elevation = cor.test(catchment_var_dt$Overall_RatioDownUp_NDVI, catchment_var_dt$Elevation, method = 'pearson')$estimate, 
  Slope = cor.test(catchment_var_dt$Overall_RatioDownUp_NDVI, catchment_var_dt$Slope, method = 'pearson')$estimate, 
  Aspect = cor.test(catchment_var_dt$Overall_RatioDownUp_NDVI, catchment_var_dt$Aspect, method = 'pearson')$estimate, 
  Latitude = cor.test(catchment_var_dt$Overall_RatioDownUp_NDVI, catchment_var_dt$Latitude, method = 'pearson')$estimate, 
  HLI = cor.test(catchment_var_dt$Overall_RatioDownUp_NDVI, catchment_var_dt$HLI, method = 'pearson')$estimate, 
  Pct_evergreen = cor.test(catchment_var_dt$Overall_RatioDownUp_NDVI, catchment_var_dt$pct_evergreen, method = 'pearson')$estimate,
  AreaSqKm = cor.test(catchment_var_dt$Overall_RatioDownUp_NDVI, catchment_var_dt$AreaSqKm_x, method = 'pearson')$estimate 
)

ratio_row <- data.table(
  Predictor = 'Downslope:Upslope NDVI Ratio', 
  Tmin = cor.test(catchment_var_dt$Overall_RatioDownUp_NDVI, catchment_var_dt$tmin, method = 'pearson')$p.value, 
  Tmax = cor.test(catchment_var_dt$Overall_RatioDownUp_NDVI, catchment_var_dt$tmax, method = 'pearson')$p.value, 
  VPD = cor.test(catchment_var_dt$Overall_RatioDownUp_NDVI, catchment_var_dt$vpd, method = 'pearson')$p.value, 
  Precipitation = cor.test(catchment_var_dt$Overall_RatioDownUp_NDVI, catchment_var_dt$prcp, method = 'pearson')$p.value, 
  PAR = cor.test(catchment_var_dt$Overall_RatioDownUp_NDVI, catchment_var_dt$par, method = 'pearson')$p.value,
  Elevation = cor.test(catchment_var_dt$Overall_RatioDownUp_NDVI, catchment_var_dt$Elevation, method = 'pearson')$p.value, 
  Slope = cor.test(catchment_var_dt$Overall_RatioDownUp_NDVI, catchment_var_dt$Slope, method = 'pearson')$p.value, 
  Aspect = cor.test(catchment_var_dt$Overall_RatioDownUp_NDVI, catchment_var_dt$Aspect, method = 'pearson')$p.value, 
  Latitude = cor.test(catchment_var_dt$Overall_RatioDownUp_NDVI, catchment_var_dt$Latitude, method = 'pearson')$p.value, 
  HLI = cor.test(catchment_var_dt$Overall_RatioDownUp_NDVI, catchment_var_dt$HLI, method = 'pearson')$p.value, 
  Pct_evergreen = cor.test(catchment_var_dt$Overall_RatioDownUp_NDVI, catchment_var_dt$pct_evergreen, method = 'pearson')$p.value,
  AreaSqKm = cor.test(catchment_var_dt$Overall_RatioDownUp_NDVI, catchment_var_dt$AreaSqKm_x, method = 'pearson')$p.value 
)

ratio_row[,2:ncol(ratio_row)] <- round(ratio_row[,2:ncol(ratio_row)], 3)

trend_row <- data.table(
  Predictor = 'Trend in Downslope:Upslope NDVI Ratio', 
  Tmin = cor.test(catchment_var_dt$Trend, catchment_var_dt$tmin, method = 'pearson')$estimate, 
  Tmax = cor.test(catchment_var_dt$Trend, catchment_var_dt$tmax, method = 'pearson')$estimate, 
  VPD = cor.test(catchment_var_dt$Trend, catchment_var_dt$vpd, method = 'pearson')$estimate, 
  Precipitation = cor.test(catchment_var_dt$Trend, catchment_var_dt$prcp, method = 'pearson')$estimate, 
  PAR = cor.test(catchment_var_dt$Trend, catchment_var_dt$par, method = 'pearson')$estimate,
  Elevation = cor.test(catchment_var_dt$Trend, catchment_var_dt$Elevation, method = 'pearson')$estimate, 
  Slope = cor.test(catchment_var_dt$Trend, catchment_var_dt$Slope, method = 'pearson')$estimate, 
  Aspect = cor.test(catchment_var_dt$Trend, catchment_var_dt$Aspect, method = 'pearson')$estimate, 
  Latitude = cor.test(catchment_var_dt$Trend, catchment_var_dt$Latitude, method = 'pearson')$estimate, 
  HLI = cor.test(catchment_var_dt$Trend, catchment_var_dt$HLI, method = 'pearson')$estimate, 
  Pct_evergreen = cor.test(catchment_var_dt$Trend, catchment_var_dt$pct_evergreen, method = 'pearson')$estimate,
  AreaSqKm = cor.test(catchment_var_dt$Trend, catchment_var_dt$AreaSqKm_x, method = 'pearson')$estimate 
)

trend_row <- data.table(
  Predictor = 'Trend in Downslope:Upslope NDVI Ratio', 
  Tmin = cor.test(catchment_var_dt$Trend, catchment_var_dt$tmin, method = 'pearson')$p.value, 
  Tmax = cor.test(catchment_var_dt$Trend, catchment_var_dt$tmax, method = 'pearson')$p.value, 
  VPD = cor.test(catchment_var_dt$Trend, catchment_var_dt$vpd, method = 'pearson')$p.value, 
  Precipitation = cor.test(catchment_var_dt$Trend, catchment_var_dt$prcp, method = 'pearson')$p.value, 
  PAR = cor.test(catchment_var_dt$Trend, catchment_var_dt$par, method = 'pearson')$p.value,
  Elevation = cor.test(catchment_var_dt$Trend, catchment_var_dt$Elevation, method = 'pearson')$p.value, 
  Slope = cor.test(catchment_var_dt$Trend, catchment_var_dt$Slope, method = 'pearson')$p.value, 
  Aspect = cor.test(catchment_var_dt$Trend, catchment_var_dt$Aspect, method = 'pearson')$p.value, 
  Latitude = cor.test(catchment_var_dt$Trend, catchment_var_dt$Latitude, method = 'pearson')$p.value, 
  HLI = cor.test(catchment_var_dt$Trend, catchment_var_dt$HLI, method = 'pearson')$p.value, 
  Pct_evergreen = cor.test(catchment_var_dt$Trend, catchment_var_dt$pct_evergreen, method = 'pearson')$p.value,
  AreaSqKm = cor.test(catchment_var_dt$Trend, catchment_var_dt$AreaSqKm_x, method = 'pearson')$p.value 
)

trend_row[,2:ncol(trend_row)] <- round(trend_row[,2:ncol(trend_row)], 3)

ndvi_row <- data.table(
  Predictor = 'Mean_NDVI', 
  Tmin = cor.test(catchment_var_dt$Mean_NDVI, catchment_var_dt$tmin, method = 'pearson')$estimate, 
  Tmax = cor.test(catchment_var_dt$Mean_NDVI, catchment_var_dt$tmax, method = 'pearson')$estimate, 
  VPD = cor.test(catchment_var_dt$Mean_NDVI, catchment_var_dt$vpd, method = 'pearson')$estimate, 
  Precipitation = cor.test(catchment_var_dt$Mean_NDVI, catchment_var_dt$prcp, method = 'pearson')$estimate, 
  PAR = cor.test(catchment_var_dt$Mean_NDVI, catchment_var_dt$par, method = 'pearson')$estimate,
  Elevation = cor.test(catchment_var_dt$Mean_NDVI, catchment_var_dt$Elevation, method = 'pearson')$estimate, 
  Slope = cor.test(catchment_var_dt$Mean_NDVI, catchment_var_dt$Slope, method = 'pearson')$estimate, 
  Aspect = cor.test(catchment_var_dt$Mean_NDVI, catchment_var_dt$Aspect, method = 'pearson')$estimate, 
  Latitude = cor.test(catchment_var_dt$Mean_NDVI, catchment_var_dt$Latitude, method = 'pearson')$estimate, 
  HLI = cor.test(catchment_var_dt$Mean_NDVI, catchment_var_dt$HLI, method = 'pearson')$estimate, 
  Pct_evergreen = cor.test(catchment_var_dt$Mean_NDVI, catchment_var_dt$pct_evergreen, method = 'pearson')$estimate,
  AreaSqKm = cor.test(catchment_var_dt$Mean_NDVI, catchment_var_dt$AreaSqKm_x, method = 'pearson')$estimate 
)
ndvi_row[,2:ncol(ndvi_row)] <- round(ndvi_row[,2:ncol(ndvi_row)], 3)

fwrite(rbind(ndvi_row, ratio_row, trend_row), paste0(home, "/Data/table_ndvi_cors.csv"))


#################################################################################################################################################
#################################################################################################################################################
#################################################################################################################################################
# Make a figure for the trend in each climate variable 
# A panel with a map for each variable (tmin, tmax, prcp, vp, vpd, par)
# A final panel with the % of catchments with a significant negative vs positive trend 
catch <- st_read(paste0(home, "/Data/Catchments/Headwater/headwater_catchments_perm_forest_32617.shp"))
catch <- catch[,c("NHDPlusID", "geometry")]
roi <- st_read("/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling/Data/ROI/blue_ridge_plus_reference.shp")
roi <- st_transform(roi, crs="EPSG:32617")

# find the scale limits 
vpd <- fread(paste0(home, "/Data/Climate/Summary/gs_trends_gs_", "vpd", ".csv"))
vp <- fread(paste0(home, "/Data/Climate/Summary/gs_trends_gs_", "vp", ".csv"))
tmin <- fread(paste0(home, "/Data/Climate/Summary/gs_trends_gs_", "tmin", ".csv"))
tmax <- fread(paste0(home, "/Data/Climate/Summary/gs_trends_gs_", "tmax", ".csv"))
prcp <- fread(paste0(home, "/Data/Climate/Summary/gs_trends_gs_", "prcp", ".csv"))
par <- fread(paste0(home, "/Data/Climate/Summary/gs_trends_gs_", "par", ".csv"))

all_slopes <- c(vpd$slope, vp$slope, tmin$slope, tmax$slope, prcp$slope, par$slope)
upper_limit <- quantile(all_slopes, 0.975)
lower_limit <- quantile(all_slopes, 0.025)


plot_climate_trend <- function(VAR, title){
  
  df <- fread(paste0(home, "/Data/Climate/Summary/gs_trends_gs_", VAR, ".csv"))
  df_shp <- merge(catch, df[df$pvalue <= 0.05,], inner=TRUE, by="NHDPlusID")
  df_shp <- cbind(df_shp, st_coordinates(st_centroid(df_shp)))
  df_shp[df_shp$slope < lower_limit, 'slope'] <- lower_limit
  df_shp[df_shp$slope > upper_limit, 'slope'] <- upper_limit
  
  
  clim_plot <- ggplot() + 
    geom_sf(data=roi, fill='white', color='black')+ 
    geom_point(data=df_shp, aes(x=X, y=Y, color=slope), alpha=0.4, size=0.001) + 
    theme_classic() + 
    xlab("") + 
    ylab("") + 
    ggtitle(title) +
    #scale_color_gradient2(low = "#F62C00", mid = "light yellow", high = "#440154", midpoint = 0, name=VAR, limits=c(lower_limit, upper_limit)) + 
    #theme(legend.title = element_blank()) + 
    #theme(legend.position = "none")
    theme(legend.position = c(0.95, 0.35)) + 
    scale_color_gradient2(low = "#F62C00", mid = "light yellow", high = "#440154", midpoint = 0, name="", limits=c(lower_limit, upper_limit)) 
  
  return(clim_plot)
}

df <- fread(paste0(home, "/Data/Climate/Summary/gs_trends_gs_", "tmin", ".csv"))
df_shp <- merge(catch, df[df$pvalue <= 0.05,], inner=TRUE, by="NHDPlusID")
df_shp <- cbind(df_shp, st_coordinates(st_centroid(df_shp)))
df_shp[df_shp$slope < lower_limit, 'slope'] <- lower_limit
df_shp[df_shp$slope > upper_limit, 'slope'] <- upper_limit

legend <- ggplot() + 
  geom_sf(data=roi, fill='white', color='black')+ 
  geom_point(data=df_shp, aes(x=X, y=Y, color=slope), alpha=0.4, size=0.001) + 
  scale_color_gradient2(low = "#F62C00", mid = "light yellow", high = "#440154", midpoint = 0, name="Sen's Slope", limits=c(lower_limit, upper_limit)) + 
  theme_classic() + 
  xlab("") + 
  ylab("") + 
  ggtitle("") + 
  #theme(legend.title = element_blank()) + 
  theme(legend.position = "right", 
        legend.text=element_text(size=12)) +
  theme(axis.text=element_text(size=12, color = 'black'),
        axis.title=element_text(size=12,face="bold")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1)) 
legend <- ggpubr::get_legend(legend)

tmin_trend_plot <- plot_climate_trend("tmin", bquote(T[min]))
tmax_trend_plot <- plot_climate_trend("tmax", bquote(T[max]))
#vp_trend_plot <- plot_climate_trend('vp')
vpd_trend_plot <- plot_climate_trend('vpd', 'VPD')
prcp_trend_plot <- plot_climate_trend('prcp', 'Prcp')
par_trend_plot <- plot_climate_trend('par', 'PAR')


tiff(paste0(home, "/Figures/climate_trends_maps.tiff"), units = 'in', res=300, width = 8, height = 10)
#ggarrange(tmin_trend_plot, tmax_trend_plot, vp_trend_plot, vpd_trend_plot, prcp_trend_plot, par_trend_plot,
#          nrow = 3, ncol = 2, 
#          common.legend = TRUE,
#          legend.grob=legend,
#          legend = "right", 
#          labels='AUTO')
ggarrange(tmin_trend_plot, tmax_trend_plot, vpd_trend_plot, prcp_trend_plot, par_trend_plot,nrow = 3, ncol = 2, labels='AUTO')
dev.off()



#################################################################################################################################################
#################################################################################################################################################
#################################################################################################################################################
## Make a plot with different panels for raster maps of every predictor variable 
## Elevation, slope, aspect, tmin, tmax, prcp, vp, vpd, par 
## crop and mask each of these to the boundary 
## save to a folder 

# get the growing season climate tifs
clim_files <- list.files(paste0(home, "/Data/Climate/tifs"), pattern = "*_AMJJA.tif$", full.names=TRUE)

# rerpoject the roi to the same crs as the climate tifs
clim_temp<- raster(clim_files[1])
roi <- st_read("/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling/Data/ROI/blue_ridge_plus_reference.shp")
roi <- st_transform(roi, crs=crs("EPSG:32617"))

# mask climate tifs to roi
for(file in clim_files){
  r <- raster(file)
  r <- projectRaster(r, crs=CRS("EPSG:32617"))
  r_mask <- mask(crop(r, roi), roi)
  writeRaster(r_mask, paste0(home, "/Data/predictors_masked/", basename(file)), overwrite = TRUE)
}
vp <- raster("/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling//Data/predictors_masked/climatology_vp_AMJJA.tif")
vp <- vp/1000
writeRaster(vp, "/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling//Data/predictors_masked/climatology_vp_AMJJA.tif", overwrite = TRUE)

top_files <- c("/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling/Data/Topography/Aspect/aspect_10m_sbr.tif", 
               "/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling/Data/Topography/Elevation/elevation_10m_sbr_32617.tif", 
               "/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling/Data/Topography/Slope/slope_10m_sbr.tif")

for(file in top_files){
  r <- raster(file)
  r_mask <- mask(crop(r, roi), roi)
  writeRaster(r_mask, paste0(home, "/Data/predictors_masked/", basename(file)), overwrite = TRUE)
}

file <- '/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling/Data/nlcd_permanent_forest/nlcd_2016.tif'
r <- raster(file)
r <- projectRaster(r, crs=CRS("EPSG:32617"), method = 'ngb')
r_mask <- mask(crop(r, roi), roi)
writeRaster(r_mask, paste0(home, "/Data/predictors_masked/", basename(file)), overwrite = TRUE)

# evergreen is 42 so make 42 a color and everything else a color 
r_mask <- raster(paste0(home, "Data/predictors_masked/", basename(file)))
r_e <- r_mask
r_e[!r_e == 42] <- NA


tiff(paste0(home, "/Figures/map_predictors.tiff"), units = 'in', res=1000, width = 14, height = 10)
par(mfrow = c(3,3), omi=c(0.3,0.3,0.3,0.3), mar = c(2, 2, 0.1, 0.1)) # plt=c(0.1,0.9,0,0.7)
all_files <- list.files(paste0(home, "/Data/predictors_masked"), pattern = "*.tif$", full.names=TRUE)
final_files <- c(all_files[c(5,4,7)], all_files[3], all_files[2], all_files[c(8,10)], all_files[1])
titles <- c(expression(T[min]~(degree*C)), expression(T[max]~(degree*C)), "VPD (kPa)", "Prcp (mm)", expression(PAR~(mol~m^-2~d^-1)), "Elevation (m)", expression(Slope~(degree)), expression(Aspect~(degree)))
K <- 0

for(i in final_files){
  K <- K + 1
  r <- raster(i)
  plot(r, col = viridis(100), axes = FALSE, box = FALSE, legend.width=1.5)
  title(titles[K], line=-1)
}
plot(r_e, col = viridis(100)[20], axes = FALSE, box = FALSE, legend.width=1.5, colNA='white', legend=F)
plot(st_geometry(roi), col=rgb(red = 0, green = 0, blue = 0, alpha = 0), add=TRUE)
legend(x=380000, y=3850000, legend = c("Evergreen \n Forest"), fill = viridis(100)[20], box.lwd = 0,box.col = "white",bg = "white")
title(expression('Evergreen Forest'), line=-1)

dev.off()


tiff(paste0(home, "/Figures/Final_Figures/Fig2.tiff"), units = 'in', res=1000, width = 14, height = 10)
par(mfrow = c(3,3), omi=c(0.3,0.3,0.3,0.3), mar = c(2, 2, 0.1, 0.1)) # plt=c(0.1,0.9,0,0.7)
all_files <- list.files(paste0(home, "/Data/predictors_masked"), pattern = "*.tif$", full.names=TRUE)
final_files <- c(all_files[c(5,4,7)], all_files[3], all_files[2], all_files[c(8,10)], all_files[1])
titles <- c(expression(T[min]~(degree*C)), expression(T[max]~(degree*C)), "VPD (kPa)", "Prcp (mm)", expression(PAR~(mol~m^-2~d^-1)), "Elevation (m)", expression(Slope~(degree)), expression(Aspect~(degree)))
K <- 0

for(i in final_files){
  K <- K + 1
  r <- raster(i)
  plot(r, col = viridis(100), axes = FALSE, box = FALSE, legend.width=1.5)
  title(titles[K], line=-1, cex.main=2, family="ArialMT")
}
plot(r_e, col = viridis(100)[20], axes = FALSE, box = FALSE, legend.width=1.5, colNA='white', legend=F)
plot(st_geometry(roi), col=rgb(red = 0, green = 0, blue = 0, alpha = 0), add=TRUE)
legend(x=380000, y=3850000, legend = c("Evergreen \n Forest"), fill = viridis(100)[20], box.lwd = 0,box.col = "white",bg = "white")
title(expression('Evergreen Forest'), line=-1, cex.main=2, family="ArialMT")

dev.off()



par(bty="l")
x <- seq(2000, 2020, 1)
y <- c(-1.5, -1.8, -0.8, 0.3, 0.15, 0.9, -0.8, -1.3, -2.1, -0.5, 
       0.1, 0.3, -0.1, 1, 0.8, 1.1, -0.2, -1.1, 0.25, -1, -0.5)
plot(x, y, ylab="Forest water use anomaly", xlab="Date", main="",
     xlim=c(2000, 2020), ylim=c(-2.5, 2.5), type='b')
rect(2000, -4, 2002, 3,
     col= rgb(red = 0.7, green = 0, blue = 0, alpha = 0.2), 
     border= rgb(red = 0.7, green = 0, blue = 0, alpha = 0.2))
rect(2006, -4, 2009, 3,
     col= rgb(red = 0.7, green = 0, blue = 0, alpha = 0.2), 
     border= rgb(red = 0.7, green = 0, blue = 0, alpha = 0.2))
rect(2017, -4, 2016, 3,
     col= rgb(red = 0.7, green = 0, blue = 0, alpha = 0.2), 
     border= rgb(red = 0.7, green = 0, blue = 0, alpha = 0.2))
rect(2019, -4, 2020, 3,
     col= rgb(red = 0.7, green = 0, blue = 0, alpha = 0.2), 
     border= rgb(red = 0.7, green = 0, blue = 0, alpha = 0.2))



