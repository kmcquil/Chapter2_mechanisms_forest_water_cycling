#############################################################################################################################
#############################################################################################################################
## Calculate trends in low flows, ET, ET:P, and make figures of all of the flow dynamics including recession slope coefficients 

#library(devtools)
#devtools::install_github("https://github.com/cran/lfstat")
library(dataRetrieval)
library(data.table)
library(sf)
library(lubridate)
library(lfstat)
library(wql)
library(kableExtra)
library(tidyr)
library(cowplot)
library(ggpubr)
library(viridis)
library(colorspace)
library(RColorBrewer)

home <- "/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling"
roi <- st_read(paste0(home, "/Data/Catchments/Reference/gages_ii/reference_keep.shp"))

get_baseflow <- function(gage){
  # set up the data table 
  result <- data.table(Gage = gage, 
                       BFI_overall = as.numeric(), 
                       BFI_intercept = as.numeric(), 
                       BFI_slope = as.numeric(), 
                       BFI_pvalue = as.numeric(), 
                       MAM7_overall = as.numeric(), 
                       MAM7_slope = as.numeric(), 
                       MAM_pvalue = as.numeric(), 
                       BFI_summer_overall = as.numeric(), 
                       BFI_summer_intercept = as.numeric(),
                       BFI_summer_slope = as.numeric(), 
                       BFI_summer_pvalue = as.numeric(), 
                       MAM7_summer_overall = as.numeric(), 
                       MAM7_summer_intercept = as.numeric(), 
                       MAM7_summer_slope = as.numeric(), 
                       MAM7_summer_pvalue = as.numeric())
  # download streamflow data 
  discharge <- readNWISdv(siteNumbers=gage, startDate = '1984-01-01', endDate = '2009-12-31', parameterCd = '00060') # mean daily discharge in ft3/s
  #discharge <- discharge[substr(discharge$Date, 6, 7) %in% c("08","09", "10"),]
  discharge <- as.data.table(discharge)
  colnames(discharge) <- c("Agency", "Site", "Date", "flow", "QC")
  #discharge[!QC == "A", flow := NA]
  discharge <- na.omit(discharge, cols="flow")
  discharge_zoo <- zoo(discharge$flow, discharge$Date)
  
  # create the low flow object that the package uses 
  lf_obj <-as.lfobj(discharge_zoo, hyearstart=10)
  lf_obj <- createlfobj(lf_obj, baseflow = TRUE, meta = NULL)
  
  ###################################################################################################
  ## Full year 
  # calculate the baseflow index and find the trend with sens slope and mann kendall trend test 
  # overall bfi 
  result$BFI_overall[1] <-BFI(lf_obj, year = "any", breakdays = NULL, yearly = FALSE)
  # annual bfi and trend 
  bfi <- BFI(lf_obj, year = "any", breakdays = NULL, yearly = TRUE)
  bfi <- as.data.frame(bfi)
  bfi <- setDT(bfi, keep.rownames = TRUE)[]
  colnames(bfi) <- c("Year", "BFI")
  bfi$Year <- as.numeric(bfi$Year)
  
  sens_bfi <- mannKen(x=bfi$BFI)
  result$BFI_slope <- sens_bfi$sen.slope
  result$BFI_pvalue <- sens_bfi$p.value
  #result$BFI_intercept <- median(bfi$BFI) - (sens_bfi$sen.slope * median(bfi$Year))
  result$BFI_intercept <- median(bfi$BFI) - (sens_bfi$sen.slope * median(seq(0, 25, 1)))
  # https://stats.stackexchange.com/questions/50587/intercept-calculation-in-theil-sen-estimator

  # calculate the mean annual minimum - 7 day 
  # overall 
  result$MAM7_overall <- MAM(lf_obj, n = 3, year = "any",breakdays = NULL, yearly = FALSE)
  mam7 <- MAM(lf_obj, n = 3, year = "any",breakdays = NULL, yearly = TRUE)
  sens_mam <- mannKen(x=mam7$MAn)
  result$MAM7_slope <- sens_mam$sen.slope
  result$MAM_pvalue <- sens_mam$p.value
  
  ###################################################################################################
  ## Summer (June - August )
  result$BFI_summer_overall <- as.numeric(BFI(lf_obj, year = "any", breakdays = "01/06", yearly = FALSE)[1])
  bfi <- BFI(lf_obj, year = "any", breakdays = "01/06", yearly = TRUE)
  bfi <- bfi[seq(2, 54, 2)]
  bfi <- as.data.frame(bfi)
  bfi <- setDT(bfi, keep.rownames = TRUE)
  bfi[,1] <- substr(bfi$rn, 25, 30)
  bfi <- na.omit(bfi)
  colnames(bfi) <- c("Year", "BFI")
  bfi$Year <- as.numeric(bfi$Year)
  sens_bfi <- mannKen(x=bfi$BFI)
  result$BFI_summer_slope <- sens_bfi$sen.slope
  result$BFI_summer_pvalue <- sens_bfi$p.value
  result$BFI_summer_intercept <- median(bfi$BFI) - (sens_bfi$sen.slope * median(seq(0, 25, 1)))
  
  # calculate the mean annual minimum - 7 day 
  # overall 
  result$MAM7_summer_overall <- MAM(lf_obj, n = 3, year = "any",breakdays = "01/08", yearly = FALSE)[1,2]
  mam7 <- MAM(lf_obj, n = 3, year = "any",breakdays = "01/08", yearly = TRUE)
  mam7 <- mam7[seq(2, nrow(mam7), 2),2:3]
  
  area_m2 <- roi[roi$GAGE_ID == gage,]$AREA
  m2_to_ft2 <- 10.7639
  ft_to_mm <- 304.8
  s_to_day <- 60*60*24
  mam7$MAn <- mam7$MAn * ((ft_to_mm * s_to_day)/(area_m2 * m2_to_ft2))
  
  sens_mam <- mannKen(x=mam7$MAn)
  result$MAM7_summer_slope <- sens_mam$sen.slope
  result$MAM7_summer_pvalue <- sens_mam$p.value
  result$MAM7_summer_intercept <- median(mam7$MAn) - (sens_mam$sen.slope * median(seq(0, 25, 1)))
  
  return(result)
  
}

# bring in the stream flow recession results to get the list of gages for the baseflow analysis 
recession_results <- fread(paste0(home, "/Data/Streamflow/baseflow_recession.csv"), 
                           colClasses=c(gage="character",A="numeric",B="numeric", slope_A='numeric', pvalue_A='numeric', 
                                        slope_B='numeric', pvalue_B='numeric'))
for(i in 1:nrow(recession_results)){
  
  result <- get_baseflow(recession_results$gage[i])
  if(i == 1){
    results <- result
  }
  else{
    results <- rbind(results, result)
  }
}

# merge the baseflow results with the recession results 
# all.y because the recession analysis accounted for gages with no missing data 
results <- merge(results, recession_results, by.x="Gage", by.y="gage", all.y=TRUE)


#######################################################################################
## Merge the results with averaged NDVI trends and stats on how much of a watershed the significant hw make up 
ndvi_results <- fread(paste0(home, "/Data/NDVI/reference_watershed_ndvi_trend_summary.csv"), colClasses=c(Gage="character", pct_hw = "numeric", wavg_rslopesig="numeric"))
perm_forest <- fread(paste0(home, "/Data/lcmap_permanent_forest/lcmap_perm_forest_reference_watersheds.csv"), colClasses=c(GAGE_ID='character', percent_forest='numeric'))
colnames(perm_forest) <- c("Gage", "pct_perm_forest")
perm_forest$Gage <- paste0('0', perm_forest$Gage)

results <- merge(results, ndvi_results)
results <- merge(results, perm_forest)

results_sub <- results[results$pvalue_A <= 0.05,]
plot(results_sub$slope_A, results_sub$wavg_rslopesig, xlab="Slope of Intercept log(a)", ylab="Slope ratio down:upslope NDVI", main="Cor = 0.70", col = 'red', pch = 16)
cor(results_sub$slope_A, results_sub$wavg_rslopesig)
# A is getting bigger (lower stability), 

results_sub <- results[results$pvalue_B <= 0.05,]
plot(results_sub$slope_B, results_sub$wavg_rslopesig, xlab="Slope of Slope(b)", ylab="Slope ratio down:upslope NDVI", main="Cor = -0.42", col = 'red', pch = 16)
cor(results_sub$slope_B, results_sub$wavg_rslopesig)

hist(results$pct_hw, xlab = "% of reference watershed covered by headwater catchments", main = "")
hist(results$pct_perm_forest, xlab = "% of reference watershed permanently forested", main = "")


## ugh fucking around 
a_ga <- results[results$pvalue_A <= 0.05 & results$slope_A > 0,]$Gage
b_ga <- results[results$pvalue_B <= 0.05 & results$slope_B < 0,]$Gage
m_ga <- results[results$MAM7_summer_pvalue <= 0.05 & results$MAM7_summer_slope < 0,]$Gage
et_ga <- flux_trends[flux_trends$ET_ratio_pvalue <= 0.05 & flux_trends$ET_ratio_slope > 0,]$Gage
e_ga <- flux_trends[flux_trends$ET_pvalue <= 0.05 & flux_trends$ET_slope > 0,]$Gage
all_ga <- c(a_ga, b_ga, m_ga, et_ga, e_ga)
length(unique(all_ga))
mean(roi[roi$GAGE_ID %in% unique(all_ga),]$AREA)/1000000

checkitout <- as.data.table(table(all_ga))
checkitout[,.(check=.N), .(N)]

# 24 watersheds have a trend in at least one metric
# 11 watersheds have a trend in just one metric 
# 7 watersheds have a trend in two metrics 
# 4 watersheds have a trend in 3 metrics 
# 2 watersheds have a trend in four metrics 

####################################################################################
## connect streamflow trends with watershed characteristics 
ref_vars <- fread("/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling/Data/Catchments/Reference/gages_ii/summarize_ws_vars.csv")
ref_vars <- ref_vars[,2:ncol(ref_vars)]
ref_vars$GAGE_ID <- as.character(paste0("0", ref_vars$GAGE_ID))
errthang <- merge(results, ref_vars, by.x="Gage", by.y="GAGE_ID")

rwf_proj <- st_transform(roi, CRS=CRS("EPSG:32617"))
rwf_proj$area_km2 <- rwf_proj$AREA/(1000000)
rwf_proj <- st_drop_geometry(rwf_proj)
rwf_proj$GAGE_ID <- as.character(rwf_proj$GAGE_ID)
errthang <- merge(errthang, rwf_proj[, c("GAGE_ID", "area_km2")], by.x="Gage", by.y="GAGE_ID")

plot_ref_vars <- function(ref_char, var, var_pvalue, xlab, ylab){
  
  df <- errthang[get(var_pvalue) <= 0.05,]
  mod <- lm(df[,get(var)] ~ df[,get(ref_char)])
  pval <- paste0("p=", round(summary(mod)$coefficients[2,4], 4))
  print(pval)
  ggplot(data=errthang[get(var_pvalue) <= 0.05,], aes(x=get(ref_char), y=get(var))) + 
    geom_point(data=errthang, aes(x=get(ref_char), y=get(var)),shape=1)+ 
    geom_point(color='forest green',size=2) +
    geom_smooth(method = "lm", se=FALSE) +
    stat_regline_equation(label.y.npc ="top", label.x.npc ="left", aes(label = ..rr.label..)) + 
    xlab(xlab) + 
    ylab(ylab) + 
    theme_classic()
  
}

elevation_a <- plot_ref_vars("Elevation","slope_A", "pvalue_A", "Elevation (m)", bquote(Trend~log(italic(a))))  
slope_a <- plot_ref_vars("Slope","slope_A", "pvalue_A", bquote(Slope~(degree)), bquote(Trend~log(italic(a))))  
aspect_a <- plot_ref_vars("Aspect","slope_A", "pvalue_A", bquote(Aspect~(degree)), bquote(Trend~log(italic(a))))  
latitude_a <- plot_ref_vars("Latitude","slope_A", "pvalue_A", "Latitude", bquote(Trend~log(italic(a))))  
prcp_a <- plot_ref_vars("Prcp","slope_A", "pvalue_A", "Prcp (mm)", bquote(Trend~log(italic(a))))  
tmin_a <- plot_ref_vars("Tmin","slope_A", "pvalue_A", bquote(italic(T[min])~(degree*C)), bquote(Trend~log(italic(a))))  
tmax_a <- plot_ref_vars("Tmax","slope_A", "pvalue_A", bquote(italic(T[max])~(degree*C)), bquote(Trend~log(italic(a))))  
vp_a <- plot_ref_vars("VP","slope_A", "pvalue_A", "VP (kPa)", bquote(Trend~log(italic(a))))  
vpd_a <- plot_ref_vars("VPD","slope_A", "pvalue_A", "VPD (kPa)", bquote(Trend~log(italic(a))))  
par_a <- plot_ref_vars("PAR","slope_A", "pvalue_A", bquote(PAR~(mol~m^-2~d^-1)), bquote(Trend~log(italic(a))))  
area_a <- plot_ref_vars("area_km2","slope_A", "pvalue_A", bquote(Area~(km^2)), bquote(Trend~log(italic(a))))  
pctfor_a <- plot_ref_vars("pct_perm_forest","slope_A", "pvalue_A", "% Permanent Forest", bquote(Trend~log(italic(a))))  

tiff(paste0(home, "/Figures/ref_watershed_relationships_a.tiff"), units = 'in', res=300, width = 9, height =9)
ggarrange(elevation_a, slope_a, aspect_a, latitude_a, 
          prcp_a, par_a, tmin_a, tmax_a, vpd_a,area_a, 
          ncol=3, nrow=4, common.legend=FALSE)
dev.off()

elevation_b <- plot_ref_vars("Elevation","slope_B", "pvalue_B", "Elevation (m)", bquote(Trend~italic(b)))  
slope_b <- plot_ref_vars("Slope","slope_B", "pvalue_B", bquote(Slope~(degree)), bquote(Trend~italic(b)))  
aspect_b <- plot_ref_vars("Aspect","slope_B", "pvalue_B",  bquote(Aspect~(degree)), bquote(Trend~italic(b)))  
latitude_b <- plot_ref_vars("Latitude","slope_B", "pvalue_B", "Latitude", bquote(Trend~italic(b)))  
prcp_b <- plot_ref_vars("Prcp","slope_B", "pvalue_B","Prcp (mm)", bquote(Trend~italic(b)))  
tmin_b <- plot_ref_vars("Tmin","slope_B", "pvalue_B",  bquote(italic(T[min])~(degree*C)), bquote(Trend~italic(b)))  
tmax_b <- plot_ref_vars("Tmax","slope_B", "pvalue_B",  bquote(italic(T[max])~(degree*C)), bquote(Trend~italic(b)))  
vp_b <- plot_ref_vars("VP","slope_B", "pvalue_B", "VP (kPa)", bquote(Trend~italic(b)))  
vpd_b <- plot_ref_vars("VPD","slope_B", "pvalue_B", "VPD (kPa)", bquote(Trend~italic(b)))  
par_b <- plot_ref_vars("PAR","slope_B", "pvalue_B", bquote(PAR~(mol~m^-2~d^-1)), bquote(Trend~italic(b)))  
area_b <- plot_ref_vars("area_km2","slope_B", "pvalue_B",bquote(Area~(km^2)), bquote(Trend~italic(b)))  
pctfor_b <- plot_ref_vars("pct_perm_forest","slope_B", "pvalue_B", "% Permanent Forest", bquote(Trend~italic(b)))  

tiff(paste0(home, "/Figures/ref_watershed_relationships_b.tiff"), units = 'in', res=300, width = 9, height =9)
ggarrange(elevation_b, slope_b, aspect_b, latitude_b, 
          prcp_b, par_b, tmin_b, tmax_b,vpd_b, area_b, 
          ncol=3, nrow=4, common.legend=FALSE)
dev.off()

#############################################################################################
## Calculate ET using a vegetation year from June - May
## I only used reference watersheds with no missing data so this should be 

calculate_runoff <- function(gage){
  
  discharge <- readNWISdv(siteNumbers=gage, startDate = '1984-01-01', endDate = '2009-12-31', parameterCd = '00060') # mean daily discharge in ft3/s
  discharge <- as.data.table(discharge)
  colnames(discharge) <- c("Agency", "Site", "Date", "flow", "QC")
  discharge <- na.omit(discharge, cols="flow")
  discharge$Date <- as.Date(discharge$Date)
  discharge$WaterYear <- ifelse(month(discharge$Date) < 6, year(discharge$Date), year(discharge$Date) + 1)
  discharge <- discharge[!WaterYear %in% c(2022, 1984), ]
  discharge <- discharge[!WaterYear > 2009,]
  
  area_m2 = as.numeric(st_drop_geometry(roi[roi$GAGE_ID == gage, 'AREA']))
  m2_to_ft2 <- 10.7639
  ft_to_mm <- 304.8
  s_to_day <- 60*60*24
  convert <- ((ft_to_mm * s_to_day) / (area_m2 * m2_to_ft2))
  discharge$runoff <- discharge$flow*convert
  
  wy_runoff <- discharge[,.(Runoff=sum(runoff)), .(Site, WaterYear)]
  colnames(wy_runoff) <- c("GAGE_ID", "WaterYear", "Runoff")
  return(wy_runoff)
  
}

# calculate the trend in runoff ratio and the trend in ET for each catchment 
et_trends <- function(gage){
  sub<- wy_flux[wy_flux$GAGE_ID == gage,]
  subc <- sub
  setorder(subc, cols='WaterYear')
  
  result <- data.table(Gage=gage, 
                       Avg_ET = mean(subc$ET), 
                       Avg_ET_ratio = mean(subc$ET_Ratio),
                       Avg_Runoff = mean(subc$Runoff), 
                       Avg_Runoff_Ratio = mean(subc$Runoff_Ratio), 
                       
                       ET_intercept = as.numeric(), 
                       ET_slope = as.numeric(), 
                       ET_pvalue = as.numeric(), 
                       
                       ET_ratio_intercept = as.numeric(), 
                       ET_ratio_slope = as.numeric(), 
                       ET_ratio_pvalue = as.numeric(), 
                       
                       Runoff_intercept = as.numeric(), 
                       Runoff_slope = as.numeric(), 
                       Runoff_pvalue = as.numeric(), 
                       
                       Runoff_ratio_intercept = as.numeric(), 
                       Runoff_ratio_slope = as.numeric(), 
                       Runoff_ratio_pvalue = as.numeric(), 
                       
                       Avg_Precip = mean(subc$Precipitation),
                       Precip_intercept = as.numeric(), 
                       Precip_slope = as.numeric(), 
                       Precip_pvalue = as.numeric(), 
                       
                       Avg_ETRunoff = mean(subc$ET/subc$Runoff), 
                       ETRunoff_intercept = as.numeric(), 
                       ETRunoff_slope = as.numeric(), 
                       ETRunoff_pvalue = as.numeric()
  )
  
  sens_et <- mannKen(x=subc$ET)
  result$ET_intercept <- median(subc$ET) - (sens_et$sen.slope * median(seq(0, 36, 1)))
  result$ET_slope <- sens_et$sen.slope
  result$ET_pvalue <- sens_et$p.value
  
  sens_et <- mannKen(x=subc$ET_Ratio)
  result$ET_ratio_intercept <- median(subc$ET_Ratio) - (sens_et$sen.slope * median(seq(0, 36, 1)))
  result$ET_ratio_slope <- sens_et$sen.slope
  result$ET_ratio_pvalue <- sens_et$p.value
  
  sens_runoff <- mannKen(x=subc$Runoff)
  result$Runoff_intercept <- median(subc$Runoff) - (sens_runoff$sen.slope * median(seq(0, 36, 1)))
  result$Runoff_slope <- sens_runoff$sen.slope
  result$Runoff_pvalue <- sens_runoff$p.value
  
  sens_runoff <- mannKen(x=subc$Runoff_Ratio)
  result$Runoff_ratio_intercept <- median(subc$Runoff_Ratio) - (sens_runoff$sen.slope * median(seq(0, 36, 1)))
  result$Runoff_ratio_slope <- sens_runoff$sen.slope
  result$Runoff_ratio_pvalue <- sens_runoff$p.value
  
  sens_p <- mannKen(x=subc$Precipitation)
  result$Precip_intercept <- median(subc$Precipitation) - (sens_p$sen.slope * median(seq(0, 36, 1)))
  result$Precip_slope <- sens_p$sen.slope
  result$Precip_pvalue <- sens_p$p.value
  
  sens_p <- mannKen(x=subc$Precipitation)
  result$Precip_intercept <- median(subc$Precipitation) - (sens_p$sen.slope * median(seq(0, 36, 1)))
  result$Precip_slope <- sens_p$sen.slope
  result$Precip_pvalue <- sens_p$p.value
  
  sens_p <- mannKen(x=(subc$ET/subc$Runoff))
  result$ETRunoff_intercept <- median(subc$ET/subc$Runoff) - (sens_p$sen.slope*median(seq(0, 36, 1)))
  result$ETRunoff_slope <- sens_p$sen.slope
  result$ETRunoff_pvalue <- sens_p$p.value
  
  return(result)
  
}

# Calcualte the runoff 
for(i in 1:nrow(recession_results)){
  runoff <- calculate_runoff(recession_results$gage[i])
  if(i ==1){
    wy_runoff <- runoff
  }else{
    wy_runoff <- rbind(wy_runoff, runoff)
  }
}
wy_runoff$WaterYear <- as.numeric(wy_runoff$WaterYear)

# Calculate the precipitation by taking the average from daymet data 
wy_precip <- fread(paste0(home, "/Data/Climate/Summary/ref_watersheds_veg_wy_precipitation.csv"))
wy_precip <- wy_precip[,2:ncol(wy_precip)]
wy_precip$GAGE_ID <- paste0(0, wy_precip$GAGE_ID)
wy_precip <- as.data.table(gather(wy_precip, WaterYear, Precipitation, `1985`:`2021`))
wy_precip$WaterYear <- as.numeric(wy_precip$WaterYear)
wy_precip <- wy_precip[wy_precip$GAGE_ID %in% recession_results$gage,]

# Calculate the ET 
wy_flux <- merge(wy_runoff, wy_precip, by=c('GAGE_ID', 'WaterYear'))
wy_flux$ET <- wy_flux$Precipitation - wy_flux$Runoff
wy_flux$Runoff_Ratio <- wy_flux$Runoff/wy_flux$Precipitation
wy_flux$ET_Ratio <- wy_flux$ET/wy_flux$Precipitation

# Calculate the trends in ET 
all_gages <- unique(wy_flux$GAGE_ID)
K = 0
for(gage in all_gages){
  K = K + 1
  if(K == 1){
    flux_trends <- et_trends(gage)
  }else{
    flux_trends <- rbind(flux_trends, et_trends(gage))
  }
}

flux_trends <- merge(flux_trends, results[,c("Gage", "pct_perm_forest")], by = "Gage")


#####################################################################################
## Map watersheds with trends 
sbr  <- st_read("/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling/Data/ROI/blue_ridge_plus_reference.shp")
sbr <- st_transform(sbr, CRS=crs("EPSG:4326"))
rwf <- roi[roi$pct_stream == 100,]
rwf <- st_transform(rwf, CRS=crs("EPSG:4326"))

## make nicer ggplots 
mycol <- rgb(0, 0, 255, max = 255, alpha = 0, names = "blue50")
rwf_a <- merge(rwf, results[results$pvalue_A <= 0.05, c("slope_A", "Gage")], inner=TRUE, by.x='GAGE_ID', by.y='Gage')
slope_a <- ggplot() + 
  geom_sf(data=sbr, fill="white", color="black", size=0.4) + 
  geom_sf(data=rwf, fill = "light grey", color='black', size=0.2) + 
  geom_sf(data=rwf_a, aes(fill=slope_A), color='black', size = 0.2)+
  scale_fill_viridis_c(name=expression(paste("Intercept log(", italic("a"), ")"))) +
  theme_void() + 
  theme(legend.position = c(0.75, 0.2)) 


rwf_b <- merge(rwf, results[results$pvalue_B <= 0.05, c("slope_B", "Gage")], inner=TRUE, by.x='GAGE_ID', by.y='Gage')
slope_b <- ggplot() + 
  geom_sf(data=sbr, fill="white", color="black", size=0.4) + 
  geom_sf(data=rwf, fill = "light grey", color='black', size=0.2) + 
  geom_sf(data=rwf_b, aes(fill=slope_B), color='black', size = 0.2)+
  scale_fill_viridis_c(name=expression(paste("Slope (", italic("b"), ")"))) +
  theme_void() + 
  theme(legend.position = c(0.75, 0.2)) 

rwf_bfi <- merge(rwf, results[results$BFI_summer_pvalue <= 0.05, c("BFI_summer_slope", "Gage")], inner=TRUE, by.x='GAGE_ID', by.y='Gage')
slope_bfi <- ggplot() + 
  geom_sf(data=sbr, fill="white", color="black", size=0.4) + 
  geom_sf(data=rwf, fill = "light grey", color='black', size = 0.2) + 
  geom_sf(data=rwf_bfi, aes(fill=BFI_summer_slope), color='black', size = 0.2)+
  scale_fill_viridis_c(name = "Trend BFI") +
  theme_void() + 
  theme(legend.position = c(0.75, 0.2)) 

rwf_mam7 <- merge(rwf, results[results$MAM7_summer_pvalue <= 0.05, c("MAM7_summer_slope", "Gage")], inner=TRUE, by.x='GAGE_ID', by.y='Gage')
slope_mam7 <- ggplot() + 
  geom_sf(data=sbr, fill="white", color="black", size=0.4) + 
  geom_sf(data=rwf, fill = "light grey", color='black', size = 0.2) + 
  geom_sf(data=rwf_mam7, aes(fill=MAM7_summer_slope), color='black', size = 0.2)+
  scale_fill_viridis_c(name = bquote(MA[7]~(mm~d^-1))) +
  theme_void() + 
  theme(legend.position = c(0.75, 0.2)) 

rwf_et <- merge(rwf, flux_trends[flux_trends$ET_pvalue <=0.05, c('Gage', 'ET_slope')], inner=TRUE, by.x='GAGE_ID', by.y='Gage')
slope_et <- ggplot() + 
  geom_sf(data=sbr, fill="white", color="black", size=0.4) + 
  geom_sf(data=rwf, fill = "light grey", color='black', size = 0.2) + 
  geom_sf(data=rwf_et, aes(fill=ET_slope), color='black', size = 0.2)+
  scale_fill_viridis_c(name = bquote(ET~(mm~y^-1))) +
  theme_void() + 
  theme(legend.position = c(0.75, 0.2))

rwf_et_ratio <- merge(rwf, flux_trends[flux_trends$ET_ratio_pvalue <=0.05, c('Gage', 'ET_ratio_slope')], inner=TRUE, by.x='GAGE_ID', by.y='Gage')
slope_et_ratio <- ggplot() + 
  geom_sf(data=sbr, fill="white", color="black", size=0.4) + 
  geom_sf(data=rwf, fill = "light grey", color='black', size = 0.2) + 
  geom_sf(data=rwf_et_ratio, aes(fill=ET_ratio_slope), color='black', size = 0.2)+
  scale_fill_viridis_c(name = "ET:P") +
  theme_void() + 
  theme(legend.position = c(0.75, 0.2)) 


tiff(paste0(home, "/Figures/low_flow_maps.tiff"), units = 'in', res=300, width = 9, height =7)
ggarrange(slope_a, slope_b, slope_et, slope_et_ratio, slope_mam7, nrow=2, ncol=3, common.legend = FALSE)
dev.off()


########################################################################################################
# Create a figure of the trend lines for each reference watershed 

m_non_sig <- results[results$MAM7_summer_pvalue > 0.05,]
m_sig <- results[results$MAM7_summer_pvalue <= 0.05,]

a_sig <- recession_results[recession_results$pvalue_A <= 0.05,]
a_non_sig <- recession_results[recession_results$pvalue_A > 0.05,]

b_sig <- recession_results[recession_results$pvalue_B <= 0.05,]
b_non_sig <- recession_results[recession_results$pvalue_B > 0.05,]

et_sig <- flux_trends[flux_trends$ET_pvalue <= 0.05, ]
et_non_sig <- flux_trends[flux_trends$ET_pvalue > 0.05, ]

et_ratio_sig <- flux_trends[flux_trends$ET_ratio_pvalue <= 0.05, ]
et_ratio_non_sig <- flux_trends[flux_trends$ET_ratio_pvalue > 0.05, ]

num_unique_gages <- length(unique(c(a_sig$gage, b_sig$gage, et_sig$Gage, et_ratio_sig$Gage, m_sig$Gage)))
unique_gages <- unique(c(a_sig$gage, b_sig$gage,et_sig$Gage, et_ratio_sig$Gage, m_sig$Gage))
q4 <- qualitative_hcl(num_unique_gages, palette = "Dark 3")
q4 <- qualitative_hcl(num_unique_gages, palette = "Harmonic")

q4 <- colorRampPalette(brewer.pal(8, "Dark2"))(num_unique_gages)
q4 <- colorRampPalette(brewer.pal(8, "Paired"))(num_unique_gages)
color_dt <- setDT(as.list(q4))[]
colnames(color_dt) <- unique_gages

# also do line types 
#l4 <- c('02069700'="solid", '02143000'="dashed", '02178400'="dotted", '03479000'="solid", '03498500'="dashed", '02053800'="dotted", '02070000'="solid", 
#        '02152100'="dashed", '03460000'="dotted", '03463300'="solid", '03471500'="dashed", '03473000'="dotted", '03500000'="solid", '02143040'="dashed", 
#        '02177000'="dotted", '03504000'="solid")

mam7_plot <- ggplot() + 
  scale_y_continuous(limits=c(0, 1.3))+ 
  scale_x_continuous(limits=c(0, 25),  breaks = seq(0, 25, 5),labels=seq(1984, 2009, 5)) + 
  geom_abline(data=m_non_sig, aes(slope=MAM7_summer_slope, intercept=MAM7_summer_intercept),color='grey', linetype='solid',alpha = 0.5) + 
  geom_abline(data=m_sig, aes(slope=MAM7_summer_slope, intercept=MAM7_summer_intercept, color=Gage), size=1, alpha=0.7) + 
  scale_color_manual(values = color_dt, name='USGS Gage') + 
  theme_classic() + 
  theme(legend.key.size =  unit(0.3, "in")) + 
  xlab("Year") + 
  ylab(bquote(MA[7]~(mm~d^-1)))+ 
  theme(legend.position='none')


a_plot <- ggplot() + 
  scale_y_continuous(limits=c(-1.9, -0.5))+ 
  scale_x_continuous(limits=c(0, 25),  breaks = seq(0, 25, 5),labels=seq(1984, 2009, 5)) + 
  theme_classic() + 
  geom_abline(data=a_non_sig, aes(slope=slope_A, intercept=int_A),color='grey', alpha = 0.5) + 
  scale_color_manual(values = color_dt, name='USGS Gage') + 
  #scale_linetype_manual(values=l4, name='USGS Gage') + 
  #geom_abline(data=a_sig, aes(slope=slope_A, intercept=int_A, color = gage, linetype=gage), size=1.5) + 
  geom_abline(data=a_sig, aes(slope=slope_A, intercept=int_A, color = gage), size=1, alpha=0.7) + 
  xlab("Year") + 
  ylab(expression(paste("Intercept log(", italic("a"), ")")))+ 
  theme(legend.position='none')

b_plot <- ggplot() + 
  scale_y_continuous(limits=c(0.8, 3.7))+ 
  scale_x_continuous(limits=c(0, 25),  breaks = seq(0, 25, 5),labels=seq(1984, 2009, 5)) + 
  theme_classic() + 
  geom_abline(data=b_non_sig, aes(slope=slope_B, intercept=int_B),color='grey', alpha = 0.5) + 
  scale_color_manual(values = color_dt, name='USGS Gage') + 
  #scale_linetype_manual(values=l4, name='USGS Gage') + 
  #geom_abline(data= b_sig, aes(slope=slope_B, intercept=int_B, color = gage, linetype=gage), size=1.5) + 
  geom_abline(data= b_sig, aes(slope=slope_B, intercept=int_B, color = gage), size=1, alpha=0.7) + 
  xlab("Year") + 
  ylab(expression(paste("Slope (", italic("b"), ")")))+ 
  theme(legend.position='none')

et_plot <- ggplot() + 
  scale_y_continuous(limits=c(420, 1150))+ 
  scale_x_continuous(limits=c(0, 25),  breaks = seq(0, 25, 5),labels=seq(1984, 2009, 5)) + 
  theme_classic() + 
  geom_abline(data=et_non_sig, aes(slope=ET_slope, intercept=ET_intercept),color='grey', alpha = 0.5) + 
  geom_abline(data= et_sig, aes(slope=ET_slope, intercept=ET_intercept, color = Gage), size=1, alpha=0.7) + 
  scale_color_manual(values = color_dt, name='USGS Gage') + 
  #scale_linetype_manual(values=l4, name='USGS Gage') + 
  xlab("Year") + 
  ylab(bquote(ET~(mm~y^-1)))+ 
  theme(legend.position='none')

et_ratio_plot <- ggplot() + 
  scale_y_continuous(limits=c(0.3, 0.8)) + 
  scale_x_continuous(limits=c(0, 25),  breaks = seq(0, 25, 5),labels=seq(1984, 2009, 5)) + 
  theme_classic() + 
  geom_abline(data=et_ratio_non_sig, aes(slope=ET_ratio_slope, intercept=ET_ratio_intercept),color='grey', alpha = 0.5) + 
  geom_abline(data= et_ratio_sig, aes(slope=ET_ratio_slope, intercept=ET_ratio_intercept, color = Gage), size=1, alpha=0.7) + 
  scale_color_manual(values = color_dt, name='USGS Gage') + 
  #scale_linetype_manual(values=l4, name='USGS Gage') + 
  xlab("Year") + 
  ylab("ET:P") + 
  theme(legend.position='none')

legend <- ggplot() + 
  scale_y_continuous(limits=c(0.3, 0.7)) + 
  scale_x_continuous(limits=c(0, 25),  breaks = seq(0, 25, 5),labels=seq(1984, 2009, 5)) + 
  theme_classic() + 
  geom_abline(data=et_ratio_non_sig, aes(slope=ET_ratio_slope, intercept=ET_ratio_intercept),color='grey', alpha = 0.5) + 
  geom_abline(data= et_ratio_sig, aes(slope=ET_ratio_slope, intercept=ET_ratio_intercept, color = Gage), size=1) + 
  scale_color_manual(values = color_dt, name='USGS Gage') + 
  #scale_linetype_manual(values=l4, name='USGS Gage') + 
  xlab("Year") + 
  ylab("ET:P") 
leg <- get_legend(legend)


tiff(paste0(home, "/Figures/low_flow_trends.tiff"), units = 'in', res=300, width = 9, height = 7 )
ggarrange(a_plot, b_plot, mam7_plot, et_plot, et_ratio_plot,  leg,
          ncol=3, nrow=2)
dev.off()



######################################################################################################################
## Make another figure option of the streamflow trends 
## Keep the original time series in grey and then put the trend lines over them in color 
# ET figure 
# a grey line for each watershed time series 
# a colored trend line for each significant trend

# start by getting the MAM7 annual time series and standardize by area 

get_mam7_ts <- function(gage){
  # download streamflow data 
  discharge <- readNWISdv(siteNumbers=gage, startDate = '1984-01-01', endDate = '2009-12-31', parameterCd = '00060') # mean daily discharge in ft3/s
  #discharge <- discharge[substr(discharge$Date, 6, 7) %in% c("08","09", "10"),]
  discharge <- as.data.table(discharge)
  colnames(discharge) <- c("Agency", "Site", "Date", "flow", "QC")
  #discharge[!QC == "A", flow := NA]
  discharge <- na.omit(discharge, cols="flow")
  discharge_zoo <- zoo(discharge$flow, discharge$Date)
  
  # create the low flow object that the package uses 
  lf_obj <-as.lfobj(discharge_zoo, hyearstart=10)
  lf_obj <- createlfobj(lf_obj, baseflow = TRUE, meta = NULL)
  
  mam7 <- MAM(lf_obj, n = 3, year = "any",breakdays = "01/06", yearly = TRUE)
  mam7 <- mam7[seq(2, nrow(mam7), 2),2:3]
  colnames(mam7) <- c("Year", "MAM7")
  result <- mam7
  result$Gage <- rep(gage, nrow(result))
  return(result)
}

for(i in 1:nrow(recession_results)){
  
  s_dat <- get_mam7_ts(recession_results$gage[i])
  if(i == 1){
    datdat <- s_dat
  }
  else{
    datdat <- rbind(datdat, s_dat)
  }
}
datdat$Year_0 <- datdat$Year-1984
datdat$Gage <- as.factor(datdat$Gage)
datdat <- merge(datdat, st_drop_geometry(roi[, c("AREA", "GAGE_ID")]), by.x='Gage', by.y='GAGE_ID')
m2_to_ft2 <- 10.7639
ft_to_mm <- 304.8
s_to_day <- 60*60*24
datdat$MAM7a <- datdat$MAM7 * ((ft_to_mm * s_to_day)/(datdat$AREA * m2_to_ft2))

# get teh recession slope coefficients 
rsc_ts <- fread(paste0(home, '/Data/Streamflow/rolling_recession_coefficients.csv'))
rsc_ts$Gage <- paste0(0, rsc_ts$Gage)
rsc_ts$Gage <- as.factor(rsc_ts$Gage)

num_unique_gages <- length(unique(c(a_sig$gage, b_sig$gage, et_sig$Gage, et_ratio_sig$Gage, m_sig$Gage)))
unique_gages <- unique(c(a_sig$gage, b_sig$gage,et_sig$Gage, et_ratio_sig$Gage, m_sig$Gage))
q4 <- qualitative_hcl(num_unique_gages, palette = "Dark 3")
q4 <- qualitative_hcl(num_unique_gages, palette = "Harmonic")

q4 <- colorRampPalette(brewer.pal(8, "Dark2"))(num_unique_gages)
#q4 <- colorRampPalette(brewer.pal(8, "Paired"))(num_unique_gages)
color_dt <- setDT(as.list(q4))[]
colnames(color_dt) <- unique_gages


a_plot <- ggplot() + 
  theme_classic() + 
  scale_y_continuous(limits=c(min(rsc_ts$A), max(rsc_ts$A))) + 
  scale_x_continuous(limits=c(0, 25),  breaks = seq(0, 25, 5),labels=seq(1984, 2009, 5)) + 
  geom_line(data=rsc_ts, aes(x=New_Year, y=A, group=Gage), color='grey', alpha=0.5) + 
  geom_abline(data=a_sig, aes(slope=slope_A, intercept=int_A, color=gage), size=1, alpha=0.7) + 
  scale_color_manual(values = color_dt, name='gage') + 
  xlab("Year") + 
  ylab(expression(paste("Intercept log(", italic("a"), ")"))) + 
  theme(legend.position = 'none')

b_plot <- ggplot() + 
  theme_classic() + 
  scale_y_continuous(limits=c(min(rsc_ts$B), max(rsc_ts$B))) + 
  scale_x_continuous(limits=c(0, 25),  breaks = seq(0, 25, 5),labels=seq(1984, 2009, 5)) + 
  geom_line(data=rsc_ts, aes(x=New_Year, y=B, group=Gage), color='grey', alpha=0.5) + 
  geom_abline(data=a_sig, aes(slope=slope_B, intercept=int_B, color=gage), size=1, alpha=0.7) + 
  scale_color_manual(values = color_dt, name='gage') + 
  xlab("Year") + 
  ylab(expression(paste("Slope (", italic("b"), ")"))) + 
  theme(legend.position = 'none')


wy_flux$WaterYear_0 <- wy_flux$WaterYear - 1985
et_plot <- ggplot() + 
  scale_y_continuous(limits=c(min(wy_flux$ET), max(wy_flux$ET)))+ 
  scale_x_continuous(limits=c(0, 25),  breaks = seq(0, 25, 5),labels=seq(1984, 2009, 5)) + 
  theme_classic() + 
  geom_line(data=wy_flux, aes(x=WaterYear_0, y=ET, group=GAGE_ID), color='grey', alpha=0.5) + 
  geom_abline(data= et_sig, aes(slope=ET_slope, intercept=ET_intercept, color = Gage), size=1, alpha=0.7) + 
  scale_color_manual(values = color_dt, name='GAGE_ID') + 
  xlab("Year") + 
  ylab(bquote(ET~(mm~y^-1))) + 
  theme(legend.position='none')

et_ratio_plot <- ggplot() + 
  scale_y_continuous(limits=c(min(wy_flux$ET_Ratio), max(wy_flux$ET_Ratio)))+ 
  scale_x_continuous(limits=c(0, 25),  breaks = seq(0, 25, 5),labels=seq(1984, 2009, 5)) + 
  theme_classic() + 
  geom_line(data=wy_flux, aes(x=WaterYear_0, y=ET_Ratio, group=GAGE_ID), color='grey', alpha=0.5) + 
  geom_abline(data= et_ratio_sig, aes(slope=ET_ratio_slope, intercept=ET_ratio_intercept, color = Gage), size=1, alpha=0.7) + 
  scale_color_manual(values = color_dt, name='GAGE_ID') + 
  xlab("Year") + 
  ylab("ET:P") + 
  theme(legend.position='none')

mam7_plot <- ggplot() + 
  scale_y_continuous(limits=c(min(datdat$MAM7a), max(datdat$MAM7a)))+ 
  scale_x_continuous(limits=c(0, 25),  breaks = seq(0, 25, 5),labels=seq(1984, 2009, 5)) + 
  theme_classic() + 
  geom_line(data=datdat, aes(x=Year_0, y=MAM7a, group=Gage),color='grey', alpha=0.5) + 
  geom_abline(data=m_sig, aes(slope=MAM7_summer_slope, intercept=MAM7_summer_intercept, color=Gage), size=1, alpha=0.7) +   xlab("Year") + 
  scale_color_manual(values = color_dt, name='Gage') + 
  ylab(bquote(MA[7]~(mm~d^-1))) + 
  theme(legend.position='none')

ex_plot <- ggplot() + 
  geom_sf(data=sbr, fill="white", color="black", size=0.4) +
  geom_sf(data=rwf, fill='white', color='black', size=0.2) + 
  geom_sf(data=rwf[rwf$GAGE_ID %in% unique_gages,], aes(fill=GAGE_ID), color='black', size=0.2) + 
  scale_fill_manual(values=color_dt, name="Gage") + 
  theme_classic() +  
  theme(legend.position='none')
  #theme(legend.position = c(0.85, 0.4), 
  #      legend.key.height = unit(0.1, 'cm'), 
  #      legend.key.width = unit(0.1, 'cm')) + 
  #guides(fill=guide_legend(ncol=4))

legend <- ggplot() + 
  scale_y_continuous(limits=c(0.3, 0.7)) + 
  scale_x_continuous(limits=c(0, 25),  breaks = seq(0, 25, 5),labels=seq(1984, 2009, 5)) + 
  theme_classic() + 
  geom_abline(data= et_ratio_sig, aes(slope=ET_ratio_slope, intercept=ET_ratio_intercept, color = Gage), size=1) + 
  scale_color_manual(values = color_dt, name='Gage') + 
  #scale_linetype_manual(values=l4, name='USGS Gage') + 
  xlab("Year") + 
  ylab("ET:P")

leg <- get_legend(ggplot() + 
  geom_sf(data=sbr, fill="white", color="black", size=0.4) +
  geom_sf(data=rwf, fill='white', color='black', size=0.2) + 
  geom_sf(data=rwf[rwf$GAGE_ID %in% unique_gages,], aes(fill=GAGE_ID), color='black', size=0.2) + 
  scale_fill_manual(values=color_dt, name="Gage") + 
  guides(fill=guide_legend(ncol=6)) + 
  #guides(fill=guide_legend(ncol=2)) + 
  theme(legend.key.height = unit(0.4, 'cm'),
        legend.key.width = unit(0.4, 'cm')))

tiff(paste0(home, "/Figures/low_flow_trends_with_ts.tiff"), units = 'in', res=300, width = 9, height = 7 )
ggarrange(a_plot, b_plot, et_plot, et_ratio_plot, mam7_plot, leg,
          ncol=3, nrow=2)
dev.off()

tiff(paste0(home, "/Figures/low_flow_trends_with_ts_map.tiff"), units = 'in', res=1000, width = 5.5, height = 7)
ggarrange(a_plot, b_plot, et_plot, et_ratio_plot, mam7_plot, ex_plot,
          ncol=2, nrow=3, common.legend=TRUE, legend.grob=leg,legend='bottom')
dev.off()

tiff(paste0(home, "/Figures/Final_Figures/Fig7.tiff"), units = 'in', res=1000, width = 5.5, height = 7)
ggarrange(a_plot, b_plot, et_plot, et_ratio_plot, mam7_plot, ex_plot,
          ncol=2, nrow=3, common.legend=TRUE, legend.grob=leg,legend='bottom')
dev.off()


