##########################################################################################
##########################################################################################
## Linear regression of yearly ratio NDVI ~ climate data with grouping by catchment and accounting for spatial autocorrelation

library(data.table)
library(sf) # 
library(nortest) 
library(caret)
library(spdep) #
library(ape)
library(MASS)
library(gstat) # 
library(nlme)
library(lme4)
library(lmerTest)
library(MuMIn)

home <- "/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling/"

#########################################################################################
# Create climate dataframe 
# column names NHPlusID, Year, temp_son, temp_djf, temp_mam, temp_jja, temp_mamjja, temp_annual, same for each climate variable 

get_climate <- function(VAR){
  dt <- fread(paste0(home, "Data/Climate/Summary/", VAR,".csv"))
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
dryday <- get_climate("drydays")
clim_dt <- Reduce(function(x, y) merge(x, y, all=TRUE, by=c("NHDPlusID", "Year")), list(tmin, tmax, prcp, vp, dryday))

#########################################################################################
# Get the response variable annual ratio downslope:upslope NDVI  
# NHDPlusID, Year, ratio
library(data.table)
home <- "/share/klmarti3/kmcquil/Chapter2_mechanisms_forest_water_cycling/"
ndvi_files <- list.files(paste0(home, "/Data/NDVI/catchment_ratio_ndvi_results"), full.names=TRUE, pattern="*.csv")

for(i in 1:length(ndvi_files)){
  if(file.size(ndvi_files[i]) == 0L){
    next
  }
  file <- fread(ndvi_files[i])
  file$year <- as.numeric(substr(file$Date, 1, 4))
  summary <- file[,.(RatioDownUp_NDVI=mean(RatioDownUp_NDVI), Mean_NDVI = mean(Mean_NDVI), WSID=mean(WSID)), .(year)]
  if(i == 1){
    result <- summary
  }else{
    result <- rbind(result, summary)
  }
}
fwrite(result, paste0(home, "/Data/NDVI/annual_ratio_ndvi.csv"))

# merge with the climate variables
ndvi_dt <- fread(paste0(home, "/Data/NDVI/annual_ratio_ndvi.csv"))
colnames(ndvi_dt) <- c("Year", "RatioDownUp_NDVI", "Mean_NDVI","NHDPlusID")
ndvi_dt$NHDPlusID <- as.numeric(ndvi_dt$NHDPlusID)
full_dt <- merge(clim_dt, ndvi_dt, by = c("NHDPlusID", "Year"), inner = TRUE)

# get the centroid of each polygon and join by NHDPlusID 
roi <- st_read(paste0(home, "/Data/Catchments/Headwater/headwater_catchments_perm_forest_32617.shp"))
roi <- cbind(roi, st_coordinates(st_centroid(roi)))
roi <- st_drop_geometry(roi[,c("NHDPlusID", "X", "Y")])
full_dt <- merge(full_dt, roi, by = "NHDPlusID")

# calculate VIF of all climate variables 
lin_mod <- lm(RatioDownUp_NDVI ~ tmin_WY + tmin_JJA + tmin_MAM + tmin_DJF + tmin_SON + tmin_MAMJJA +  
                tmax_WY + tmax_JJA + tmax_MAM + tmax_DJF + tmax_SON + tmax_MAMJJA + 
                prcp_WY + prcp_JJA + prcp_MAM + prcp_DJF + prcp_SON + prcp_MAMJJA + 
                vp_WY + vp_JJA + vp_MAM + vp_DJF + vp_SON + vp_MAMJJA +
                drydays_WY + drydays_JJA + drydays_MAM + drydays_DJF + drydays_SON + drydays_MAMJJA + 
                X + Y, data = full_dt)


lin_mod <- lm(RatioDownUp_NDVI ~ tmin_WY + tmin_JJA + tmin_MAM + tmin_MAMJJA +  
                tmax_WY + tmax_JJA + tmax_MAM + tmax_MAMJJA + 
                prcp_WY + prcp_JJA + prcp_MAM + prcp_MAMJJA + 
                vp_WY + vp_JJA + vp_MAM + vp_MAMJJA +
                drydays_WY + drydays_JJA + drydays_MAM + drydays_MAMJJA + 
                X + Y, data = full_dt)

lin_mod <- lm(RatioDownUp_NDVI ~  tmin_MAMJJA +  + tmax_MAMJJA +  prcp_MAMJJA +  vp_MAMJJA +drydays_MAMJJA + 
                X + Y, data = full_dt)

lin_mod <- lm(RatioDownUp_NDVI ~   tmax_MAMJJA +  prcp_MAMJJA +  vp_MAMJJA +drydays_MAMJJA + 
                X + Y, data = full_dt)

lin_mod <- lm(RatioDownUp_NDVI ~  tmin_MAMJJA  + tmax_MAMJJA +  prcp_MAMJJA + drydays_MAMJJA + 
                X + Y, data = full_dt)

car::vif(lin_mod)

# fuck 
lin_mod <- lm(Mean_NDVI ~  tmin_MAMJJA +  + tmax_MAMJJA +  prcp_MAMJJA +  vp_MAMJJA +drydays_MAMJJA + 
                X + Y, data = full_dt)
summary(lin_mod)


# use a linear mixed effects model with the catchment ID as a random effect
# use a stepwise backward selection regression 


home <- "/share/klmarti3/kmcquil/Chapter2_mechanisms_forest_water_cycling/"
final_dt_notshp <- full_dt[1:1000,]
spat_struc <- list(corSpher(form =~ X + Y, nugget = TRUE), corExp(form =~ X + Y, nugget = TRUE), corGaus(form =~ X + Y, nugget = TRUE), corRatio(form =~ X + Y, nugget = TRUE))

potential_models <- foreach(i = 1:length(spat_struc))%dopar%{
  
  library(nlme)
  library(lme4)
  library(MASS)
  library(caret)
  
  .GlobalEnv$final_dt_notshp <- final_dt_notshp
  .GlobalEnv$spat_cor <- spat_struc[[i]]
  
  model <- gls(RatioDownUp_NDVI ~  tmin_MAMJJA +  + tmax_MAMJJA +  prcp_MAMJJA +  vp_MAMJJA +drydays_MAMJJA + X + Y + (1|NHDPlusID),
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



full_dt_nona<- na.omit(full_dt)
ugh = full_dt_nona[,.(count = .N), .(NHDPlusID)]

model <- gls(RatioDownUp_NDVI ~  tmin_MAMJJA + tmax_MAMJJA +  prcp_MAMJJA +  vp_MAMJJA +drydays_MAMJJA + X + Y + (1|NHDPlusID),
             data=full_dt, 
             #correlation = spat_cor,
             method="ML") 

model <- gls(RatioDownUp_NDVI ~  tmin_MAMJJA + tmax_MAMJJA +  prcp_MAMJJA +  vp_MAMJJA +drydays_MAMJJA + X + Y ,
             data=full_dt, 
             #correlation = spat_cor,
             method="ML") 

model <- lmer(RatioDownUp_NDVI ~  tmin_MAMJJA + tmax_MAMJJA +  prcp_MAMJJA +  vp_MAMJJA +drydays_MAMJJA + (1|NHDPlusID),
             data=full_dt,REML=TRUE) 






#########################################################################################################
# Stepwise backward regression 
# Consider grouping factors HW catchment and year 
# Grouping by HW and year introduces random effects to model that covariance
# This way I do not have to explicitly model spatial or temporal autocorrelation 

# First create a model that includes all variables and both potential random effects 
# Use a linear mixed effects model from lmer package 
# Even though the response is not linear, the data is large enough we don't need to worry 

# Run stepwise backward regression 
# check the pearson correlation and VIF of all of the full model 
# check the VIF of the final model 

full_model <- lmer(Mean_NDVI ~  tmin_MAMJJA + tmax_MAMJJA +  prcp_MAMJJA +  vp_MAMJJA +drydays_MAMJJA + (1|NHDPlusID), 
                   data=full_dt)


cor.test(full_dt$RatioDownUp_NDVI, full_dt$tmin_MAMJJA)
cor.test(full_dt$RatioDownUp_NDVI, full_dt$tmax_MAMJJA)
cor.test(full_dt$RatioDownUp_NDVI, full_dt$prcp_MAMJJA)


# Marginal R2 provides the variance explained only by fixed effects and 
# Conditional R2 provides the variance explained by the entire model, i.e., both fixed effects and random effects.
r.squaredGLMM(full_model)

trying <- cbind(full_dt[,1:2],scale(full_dt[,3:33]))
trying <- cbind(trying, full_dt[,34:35])
full_model <- lm(RatioDownUp_NDVI ~  tmin_MAMJJA + tmax_MAMJJA +  prcp_MAMJJA +  vp_MAMJJA +drydays_MAMJJA, 
                 data=trying)

full_model <- lm(Mean_NDVI ~  tmin_MAMJJA + tmax_MAMJJA +  prcp_MAMJJA +  vp_MAMJJA +drydays_MAMJJA + X + Y, 
                 data=full_dt[Mean_NDVI>=0.5,])
summary(full_model)



## Look at the correlations 
sub <- full_dt[,c("Mean_NDVI", "RatioDownUp_NDVI", "tmin_MAMJJA", "tmax_MAMJJA", "prcp_MAMJJA", "vp_MAMJJA", "drydays_MAMJJA")]
M <-cor(sub)
ggcorrplot(M, 
           hc.order = TRUE, 
           type = "lower",
           lab = TRUE)




##############################################################################################################
## Update 
full_dt <- as.data.frame(full_dt)
full_dt$NHDPlusID <- as.factor(full_dt$NHDPlusID)
model_ratio <- lmer(RatioDownUp_NDVI ~ tmin_MAMJJA + tmax_MAMJJA +  prcp_MAMJJA +  vp_MAMJJA +drydays_MAMJJA + (1|NHDPlusID), 
                   data=full_dt)
#step_model_ratio <- MASS::stepAIC(model_ratio, direction ='backward', trace=TRUE, data=full_dt)
step_model_ratio <- lmerTest::step(model_ratio,data=full_dt)
model_ratio_found <- lmer(RatioDownUp_NDVI ~ tmax_MAMJJA +  prcp_MAMJJA +  vp_MAMJJA + (1|NHDPlusID), data=full_dt)
car::vif(model_ratio_found)
r.squaredGLMM(model_ratio_found)
AIC(model_ratio_found)


model_ndvi <- lmer(Mean_NDVI ~  tmin_MAMJJA + tmax_MAMJJA +  prcp_MAMJJA +  vp_MAMJJA +drydays_MAMJJA + (1|NHDPlusID), data=full_dt)
#step_model_ndvi <- MASS::stepAIC(model_ndvi, direction='backward', trace=TRUE, data = full_dt)
step_model_ndvi <- lmerTest::step(model_ndvi, data=full_dt)
car::vif(model_ndvi)
model_ndvi_found <- lmer(Mean_NDVI ~  tmax_MAMJJA +  prcp_MAMJJA +  vp_MAMJJA + (1|NHDPlusID), data=full_dt)
car::vif(model_ndvi_found)
r.squaredGLMM(model_ndvi_found)
AIC(model_ndvi_found)


fwrite(full_dt, paste0(home, "/Data/NDVI/space_time_dt.csv"))
# try with random slopes
full_dt <- fread(paste0(home,"/Data/NDVI/space_time_dt.csv"))
model_ndvi <- lmer(Mean_NDVI ~  tmin_MAMJJA + tmax_MAMJJA +  prcp_MAMJJA +  vp_MAMJJA +drydays_MAMJJA + (1|NHDPlusID), data=full_dt)
#step_model_ndvi <- MASS::stepAIC(model_ndvi, direction='backward', trace=TRUE, data = full_dt)
step_model_ndvi <- lmerTest::step(model_ndvi, data=full_dt)
car::vif(model_ndvi)
model_ndvi_found <- lmer(Mean_NDVI ~  tmax_MAMJJA +  prcp_MAMJJA +  vp_MAMJJA + (1|NHDPlusID), data=full_dt)
car::vif(model_ndvi_found)
r.squaredGLMM(model_ndvi_found)
AIC(model_ndvi_found)














##################################################################################################
# Instead of a spatial regression using all catchments I am going to calculate the correlations 
# of the annual ratio NDVI with climate variables in different seasons 

# Game this out: 
# Using Tmin: SON, DJF, MAM, JJA, MAMJJA 
# Find correlation of ratio NDVI with each of those 
# Dataframe with columns NHDPlusID, correlation,pvalue, season, variable 

# From the dataframe, create some summaries of the correlations by season and variable --- make a table and histograms 
# Then check the associations of those correlations with topography and climatology 

library(data.table)
library(sf) 
library(ggplot2)
home <- "/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling/"

#########################################################################################
# Create climate dataframe 
# column names NHPlusID, Year, temp_son, temp_djf, temp_mam, temp_jja, temp_mamjja, temp_annual, same for each climate variable 

get_climate <- function(VAR){
  dt <- fread(paste0(home, "Data/Climate/Summary/", VAR,".csv"))
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
dryday <- get_climate("drydays")
clim_dt <- Reduce(function(x, y) merge(x, y, all=TRUE, by=c("NHDPlusID", "Year")), list(tmin, tmax, prcp, vp, dryday))

#########################################################################################
# Get the response variable annual ratio downslope:upslope NDVI  
# NHDPlusID, Year, ratio

## this chunk was on the hpc 
########################################################################################
library(data.table)
home <- "/share/klmarti3/kmcquil/Chapter2_mechanisms_forest_water_cycling/"
ndvi_files <- list.files(paste0(home, "/Data/NDVI/catchment_ratio_ndvi_results"), full.names=TRUE, pattern="*.csv")

for(i in 1:length(ndvi_files)){
  if(file.size(ndvi_files[i]) == 0L){
    next
  }
  file <- fread(ndvi_files[i])
  file$year <- as.numeric(substr(file$Date, 1, 4))
  summary <- file[,.(RatioDownUp_NDVI=mean(RatioDownUp_NDVI), Mean_NDVI = mean(Mean_NDVI), WSID=mean(WSID)), .(year)]
  if(i == 1){
    result <- summary
  }else{
    result <- rbind(result, summary)
  }
}
fwrite(result, paste0(home, "/Data/NDVI/annual_ratio_ndvi.csv"))

########################################################################################
# merge with the climate variables
ndvi_dt <- fread(paste0(home, "/Data/NDVI/annual_ratio_ndvi.csv"))
colnames(ndvi_dt) <- c("Year", "RatioDownUp_NDVI", "Mean_NDVI","NHDPlusID")
ndvi_dt$NHDPlusID <- as.numeric(ndvi_dt$NHDPlusID)
full_dt <- merge(clim_dt, ndvi_dt, by = c("NHDPlusID", "Year"), inner = TRUE)

# get the centroid of each polygon and join by NHDPlusID 
roi <- st_read(paste0(home, "/Data/Catchments/Headwater/headwater_catchments_perm_forest_32617.shp"))
roi <- cbind(roi, st_coordinates(st_centroid(roi)))
roi <- st_drop_geometry(roi[,c("NHDPlusID", "X", "Y")])
full_dt <- merge(full_dt, roi, by = "NHDPlusID")

# function to calculate the correlation of ratio ndvi with a specific climate variable/aggregation and then return a df of the results
vars <- c("tmin", "tmax", "prcp", "vp", "drydays")
seasons <- c("SON", "DJF", "MAM", "JJA", "MAMJJA", "WY")
ids <- unique(full_dt$NHDPlusID)

K <- 0
for(var in vars){
  for(season in seasons){
    clim_column<- paste0(var, "_", season)
    tik <- proc.time()
    for(id in ids){
      K <- K+1
      corr <- cor.test(full_dt[NHDPlusID == id,get("RatioDownUp_NDVI")], full_dt[NHDPlusID == id,get(clim_column)], method="pearson")
      if(K == 1){
        results <- data.table('NHDPlusID'=id, 
                              'Season'=season,
                              'Variable'=var,
                              'Correlation'=corr$estimate,
                              'Pvalue'=corr$p.value)
      }else{
        row = data.table('NHDPlusID'=id, 
                         'Season'=season,
                         'Variable'=var,
                         'Correlation'=corr$estimate,
                         'Pvalue'=corr$p.value)
        results <- rbind(results, row)
      }
    }
    tok <- proc.time() - tik
    tok
  }
}
fwrite(results, paste0(home, "Data/Climate_NDVI_correlations.csv"))


## slightly different version - use r2 instead of r 
K <- 0
for(var in vars){
  for(season in seasons){
    clim_column<- paste0(var, "_", season)
    tik <- proc.time()
    for(id in ids){
      K <- K+1
      corr <- lm(full_dt[NHDPlusID == id,get("RatioDownUp_NDVI")] ~ full_dt[NHDPlusID == id,get(clim_column)])
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
    tok <- proc.time() - tik
    tok
  }
}
fwrite(results, paste0(home, "/Data/Climate_NDVI_linmod_r2.csv"))


#results <- fread(paste0(home, "Data/Climate_NDVI_correlations.csv"))
results <- fread(paste0(home, "/Data/Climate_NDVI_linmod_r2.csv"))
results$Group <- paste0(results$Variable, " ", results$Season)
results$Group <- as.factor(results$Group)
results_sig <- results[results$Pvalue <= 0.05,]
results_sig[, .(count = .N/30044), .(Group)]

tab <- results_sig[as.numeric(NHDPlusID) %in% as.numeric(resp_sig_sub$NHDPlusID), .(count = .N/30044), .(Variable, Season)]
colnames(tab) <- c("Variable", "Season", "Fraction")
tab_wide <- as.data.table(pivot_wider(tab, names_from = Season, values_from=Fraction))

kbl(tab_wide, escape=F, align="c") %>% 
  kable_classic('striped', full_width=F)


ggplot(results, aes(x=Correlation)) + 
  geom_histogram(binwidth=0.1, color='sky blue', fill='sky blue') + 
  geom_vline(xintercept = 0, color='black', linetype='dashed') + 
  theme_bw() + 
  facet_wrap(~Group)

ggplot(results_sig, aes(x=Correlation)) + 
  geom_histogram(binwidth=0.1, color='sky blue', fill='sky blue') + 
  geom_vline(xintercept = 0, color='black', linetype='dashed') + 
  theme_bw() + 
  facet_wrap(~Group)

ggplot(results_sig, aes(x=Correlation)) + 
  geom_histogram(binwidth=0.1, color='sky blue', fill='sky blue') + 
  geom_vline(xintercept = 0, color='black', linetype='dashed') + 
  theme_bw() + 
  facet_wrap(~Group)

ggplot(results_sig, aes(x=Variable, y=Correlation, fill=Season)) + 
  geom_boxplot() + 
  theme_bw()

ggplot(results, aes(x=Variable, y=Correlation, fill=Season)) + 
  geom_boxplot() + 
  theme_bw()

ggplot(results_sig, aes(x=Variable, y=Correlation, fill=Season)) + 
  geom_boxplot() + 
  theme_bw()

ggplot(results[as.numeric(NHDPlusID) %in% as.numeric(resp_sig_sub$NHDPlusID),], aes(x=Variable, y=Correlation, fill=Season)) + 
  geom_boxplot() + 
  theme_bw()

ggplot(results_sig[as.numeric(NHDPlusID) %in% as.numeric(resp_sig_sub$NHDPlusID),], aes(x=Variable, y=Correlation, fill=Season)) + 
  geom_boxplot() + 
  theme_bw()


# down/up so a larger ratio would be 0.8 / 0.5  so there would be more ripairan biomass and it is more hydrologically connected 
# a smaller ratio means there would be more biomass upslope, and it would be less hydroologically connected 

# so a positive correlation with drydays: as drydays increase so does the ratio 
# a negative correlation with tmin: as tmin increases, the ratio decreases ---- so hotter places would be less hydrologically connected 

# calculate the correlation between correlations and topogrpahy and clim and make a nice table
# row for each group and the correlation between correlation and the variabl e
catch_avg_vars <- fread("/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling/Data/catchment_avg_vars.csv")
cols<-c('NHDPlusID', 'tmin', 'tmax', 'prcp', 'vp', 'drydays', 'Aspect', 
        'Slope', 'Elevation', 'Latitude', 'HLI')
catch_avg_vars <- catch_avg_vars[,..cols]
catch_avg_vars$NHDPlusID <- as.numeric(catch_avg_vars$NHDPlusID)
big_results <- merge(results, catch_avg_vars, how='left', by='NHDPlusID')

groups <- unique(big_results$Group)
K <- 0
for(i in groups){
  K = K+1
  sub_results <- big_results[big_results$Group == i,]
  if(K==1){
    corcors <- data.table('Group'=i, 
                      'Tmin'=cor.test(sub_results$Correlation, sub_results$tmin)$estimate, 
                      'Tmax'=cor.test(sub_results$Correlation, sub_results$tmax)$estimate,
                      'Prcp'=cor.test(sub_results$Correlation, sub_results$prcp)$estimate,
                      'VP'=cor.test(sub_results$Correlation, sub_results$vp)$estimate,
                      'Drydays'=cor.test(sub_results$Correlation, sub_results$drydays)$estimate,
                      'Aspect'=cor.test(sub_results$Correlation, sub_results$Aspect)$estimate,
                      'Slope'=cor.test(sub_results$Correlation, sub_results$Slope)$estimate,
                      'Elevation'=cor.test(sub_results$Correlation, sub_results$Elevation)$estimate,
                      'Latitude'=cor.test(sub_results$Correlation, sub_results$Latitude)$estimate,
                      'HLI'=cor.test(sub_results$Correlation, sub_results$HLI)$estimate)
  }else{
    row <- data.table('Group'=i, 
                      'Tmin'=cor.test(sub_results$Correlation, sub_results$tmin)$estimate, 
                      'Tmax'=cor.test(sub_results$Correlation, sub_results$tmax)$estimate,
                      'Prcp'=cor.test(sub_results$Correlation, sub_results$prcp)$estimate,
                      'VP'=cor.test(sub_results$Correlation, sub_results$vp)$estimate,
                      'Drydays'=cor.test(sub_results$Correlation, sub_results$drydays)$estimate,
                      'Aspect'=cor.test(sub_results$Correlation, sub_results$Aspect)$estimate,
                      'Slope'=cor.test(sub_results$Correlation, sub_results$Slope)$estimate,
                      'Elevation'=cor.test(sub_results$Correlation, sub_results$Elevation)$estimate,
                      'Latitude'=cor.test(sub_results$Correlation, sub_results$Latitude)$estimate,
                      'HLI'=cor.test(sub_results$Correlation, sub_results$HLI)$estimate)
    corcors <- rbind(corcors, row)
  }
}
cols<-c('Tmin', 'Tmax', 'Prcp', 'VP', 'Drydays', 'Aspect', 
        'Slope', 'Elevation', 'Latitude', 'HLI')
corcors[,(cols):=round(.SD, 2), .SDcols=cols]

library(kableExtra)
corcors[,2:ncol(corcors)] <- lapply(corcors[,2:ncol(corcors)], function(x){
  cell_spec(x, color=spec_color(x, end=0.9))
})
corcors[,2:ncol(corcors)] <- round(corcors[,2:ncol(corcors)], 2)
kbl(corcors, escape=F, align="c") %>% 
  kable_classic('striped', full_width=F)


sub_results <- big_results[big_results$Group == "tmin WY",]
ggplot(sub_results, aes(x=vp, y=Correlation)) + 
  geom_point(color='sky blue', alpha=0.3) + 
  theme_bw() 



# since these sort of sit around 0, even though those are mostly not significant, i want to plot this and see where we are getting the 
# positive v negative correlations on any given variable 

catch <- st_read(paste0(home, "/Data/Catchments/Headwater/headwater_catchments_perm_forest_32617.shp"))
catch <- catch[,c("NHDPlusID", "geometry")]
roi <- st_read("/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling/Data/ROI/blue_ridge_plus_reference.shp")
roi <- st_transform(roi, crs="EPSG:32617")

plot_clim_trends <- function(grp){
  grp_cors <- merge(catch, results[Group == grp], inner=TRUE, by="NHDPlusID")
  grp_cors <- cbind(grp_cors, st_coordinates(st_centroid(grp_cors)))
  grp_nonsig <- grp_cors[grp_cors$Pvalue > 0.05,]
  grp_sig <- grp_cors[grp_cors$Pvalue <= 0.05,]
  
  gp <- ggplot() + 
    geom_sf(data=roi, fill='white', color='black')+ 
    geom_point(data=grp_nonsig, aes(x=X, y=Y), color='grey', alpha=0.3) + 
    geom_point(data=grp_sig, aes(x=X, y=Y, color=Correlation), alpha=0.4) + 
    #geom_point(data=grp_sig, aes(x=X, y=Y, color=Correlation, size = abs(Correlation)), alpha=0.4) + 
    scale_color_gradient2() + 
    #scale_size(range=c(0, 5)) + 
    theme_classic() + 
    ggtitle(grp)
  return(gp)
}

all_groups <- unique(results$Group)
tmin1 <- plot_clim_trends(all_groups[1])
tmin2 <- plot_clim_trends(all_groups[2])
tmin3 <- plot_clim_trends(all_groups[3])
tmin4 <- plot_clim_trends(all_groups[4])

library(cowplot)
library(ggpubr)
ggarrange(tmin1, tmin2, tmin3, tmin4,nrow = 1, ncol = 4, 
          common.legend = TRUE,
          legend = "right")




# plot scatterplots of trends vs avg vars 
resp <- fread(paste0(home, "Data/NDVI/ratio_ndvi_trend_results.csv"))
resp$wsid <- as.numeric(resp$wsid)
resp_sig <- resp[metric == "trend_ratio_ndvi" & p_value <= 0.05,] # keep only the responses that were significant 
resp_sig_sub <- resp_sig[,c("wsid","slope", "mean_ndvi")]
colnames(resp_sig_sub) <- c("NHDPlusID", "Trend", "Mean_NDVI")

scatter_df <- merge(catch_avg_vars, resp_sig_sub, how='left', by="NHDPlusID")
scatter_df_long <- gather(scatter_df, "Variable", "Value", tmin:HLI)
var <- "HLI"
ggplot(scatter_df_long[scatter_df_long$Variable == var,], aes(x=Value, y=Trend)) + 
  geom_point(color='sky blue', alpha=0.3) + 
  theme_bw() + 
  ylab('Trend ratio downslope:upslope NDVI') + 
  xlab(var)

library(kableExtra)

kbl(results_sig[, .(count = .N/30044), .(Group)], escape=F, align="c") %>% 
  kable_classic('striped', full_width=F)



#########################################################
# make a figure of the average annual ratio
library(cowplot)
library(ggpubr)
library(patchwork)

catch <- st_read(paste0(home, "/Data/Catchments/Headwater/headwater_catchments_perm_forest_32617.shp"))
catch <- catch[,c("NHDPlusID", "geometry")]
roi <- st_read("/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling/Data/ROI/blue_ridge_plus_reference.shp")
roi <- st_transform(roi, crs="EPSG:32617")

ndvi_dt <- fread(paste0(home, "/Data/NDVI/annual_ratio_ndvi.csv"))
ndvi_summary <- ndvi_dt[,.(RatioDownUp_NDVI = mean(RatioDownUp_NDVI)), .(WSID)]
ndvi_summary$WSID <- as.numeric(ndvi_summary$WSID)
catch_summary <- merge(catch, ndvi_summary, by.x="NHDPlusID", by.y="WSID", inner=TRUE)

grp_cors <- cbind(catch_summary, st_coordinates(st_centroid(catch_summary)))
grp_cors <- grp_cors[!grp_cors$RatioDownUp_NDVI > quantile(grp_cors$RatioDownUp_NDVI, 0.99, na.rm=TRUE), ] # drop anomalies 
grp_cors <- grp_cors[!grp_cors$RatioDownUp_NDVI < quantile(grp_cors$RatioDownUp_NDVI, 0.01, na.rm=TRUE), ]

grp_cors$Group <- rep(" ", nrow(grp_cors))
ratio_violin_plot <- ggplot(grp_cors, aes(x=Group,y=RatioDownUp_NDVI)) +
  geom_violin() + 
  coord_flip() + 
  theme_minimal() + 
  xlab("") + 
  ylab("") + 
  theme(axis.title.y = element_text(size = rel(0.1)))
ratio_violin_plot
ratio_map <- ggplot() + 
  geom_sf(data=roi, fill='white', color='black')+ 
  geom_point(data=grp_cors, aes(x=X, y=Y, color=RatioDownUp_NDVI), alpha=0.5, size = 0.05) + 
  #geom_point(data=grp_sig, aes(x=X, y=Y, color=Correlation, size = abs(Correlation)), alpha=0.4) + 
  scale_color_gradient2(low = "red", mid = "white", high = "blue", midpoint = 1, name=bquote(NDVI[Down:Up])) + 
  #scale_size(range=c(0, 5)) + 
  theme_classic() + 
  ggtitle('') + 
  xlab("") + 
  ylab("") + 
  theme(legend.position = c(0.9, 0.3)) + 
  inset_element(ratio_violin_plot, 0.35, 0.07, 0.7, 0.30, align_to = 'full') #left bottom right top
ratio_map



#png(paste0(home, "/Figures/ratio_ndvi_map.png"), width=1200, height=900, res=200)
#dev.off()

###################################################################################
## make figure of the trend in NDVI 
resp <- fread(paste0(home, "/Data/NDVI/ratio_ndvi_trend_results.csv"))
resp <- resp[resp$metric == "trend_ratio_ndvi",]
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
trend_violin_plot

resp_sf[resp_sf$slope > quantile(resp_sf$slope, 0.99, na.rm=TRUE), ]$slope <- quantile(resp_sf$slope, 0.99, na.rm=TRUE) # drop anomalies 
resp_sf[resp_sf$slope < quantile(resp_sf$slope, 0.01, na.rm=TRUE), ]$slope <- quantile(resp_sf$slope, 0.01, na.rm=TRUE)
trend_map <- ggplot() + 
  geom_sf(data=roi, fill='white', color='black')+ 
  geom_point(data=resp_sf, aes(x=X, y=Y, color=slope), alpha=0.5, size = 0.05) + 
  #geom_point(data=grp_sig, aes(x=X, y=Y, color=Correlation, size = abs(Correlation)), alpha=0.4) + 
  scale_color_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0, name=bquote('d' * NDVI[Down:Up])) + 
  #scale_size(range=c(0, 5)) + 
  theme_classic() + 
  ggtitle('') + 
  xlab("") + 
  ylab("") + 
  theme(legend.position = c(0.9, 0.3)) + 
  inset_element(trend_violin_plot, 0.35, 0.07, 0.7, 0.30, align_to = 'full') #left bottom right top
trend_map


#tiff(paste0(home, "/Figures/ratio_trend_ndvi_maps.tiff"), units = 'in', res=800, width = 6, height = 4)
ggarrange(ratio_map, trend_map,nrow = 1, ncol = 2, common.legend = FALSE)
#dev.off()
ggsave(paste0(home, "/Figures/ratio_trend_ndvi_maps.tiff"), units = 'in', dpi=800, width = 10, height = 6)



################################################################################
## Make plots of the pearson correlation between ratio and predictors 
library(ggcorrplot)

catch_avg_vars <- fread("/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling/Data/catchment_avg_vars.csv")
cols<-c('NHDPlusID', 'tmin', 'tmax', 'prcp', 'vp', 'drydays', 'Aspect', 
        'Slope', 'Elevation', 'Latitude', 'HLI')
catch_avg_vars <- catch_avg_vars[,..cols]
catch_avg_vars$NHDPlusID <- as.numeric(catch_avg_vars$NHDPlusID)
catch_avg_vars <- cbind(catch_avg_vars[,c("NHDPlusID")], scale(catch_avg_vars[,2:ncol(catch_avg_vars)], center = TRUE, scale = TRUE))
ratio_results <- merge(ndvi_summary, catch_avg_vars, how='left', by.x='WSID', by.y="NHDPlusID")
ratio_results$RatioDownUp_NDVI <- scale(ratio_results$RatioDownUp_NDVI, center=TRUE, scale=TRUE)

ratio_row <- data.table(
  Predictor= 'Ratio NDVI',
  Tmin=cor.test(ratio_results$RatioDownUp_NDVI, ratio_results$tmin, method = 'pearson')$estimate,
  Tmax=cor.test(ratio_results$RatioDownUp_NDVI, ratio_results$tmax, method = 'pearson')$estimate, 
  VP=cor.test(ratio_results$RatioDownUp_NDVI, ratio_results$vp, method = 'pearson')$estimate,
  Prcp=cor.test(ratio_results$RatioDownUp_NDVI, ratio_results$prcp, method = 'pearson')$estimate, 
  Elevation=cor.test(ratio_results$RatioDownUp_NDVI, ratio_results$Elevation, method = 'pearson')$estimate, 
  Slope=cor.test(ratio_results$RatioDownUp_NDVI, ratio_results$Slope, method = 'pearson')$estimate,
  Aspect=cor.test(ratio_results$RatioDownUp_NDVI, ratio_results$Aspect, method = 'pearson')$estimate,
  Latitude=cor.test(ratio_results$RatioDownUp_NDVI, ratio_results$Latitude, method = 'pearson')$estimate 
)
ratio_row[,2:ncol(ratio_row)] <- round(ratio_row[,2:ncol(ratio_row)], 3)


ggcorrplot(round(cor(ratio_results[,2:ncol(ratio_results)]), 2),
           hc.order = TRUE,
           type = "lower",
           outline.color = "white")



trend_results <- merge(resp, catch_avg_vars, how='left', by.x='wsid', by.y="NHDPlusID")
trend_results$slope <- scale(trend_results$slope, center=TRUE, scale=TRUE)
ggcorrplot(round(cor(trend_results[,2:ncol(trend_results)]), 2),
           hc.order = TRUE,
           type = "lower",
           outline.color = "white",
           lab=TRUE)

trend_row <- data.table(
  Predictor= 'Trend ratio ndvi',
  Tmin=cor.test(trend_results$slope, trend_results$tmin, method = 'pearson')$estimate,
  Tmax=cor.test(trend_results$slope, trend_results$tmax, method = 'pearson')$estimate, 
  VP=cor.test(trend_results$slope, trend_results$vp, method = 'pearson')$estimate,
  Prcp=cor.test(trend_results$slope, trend_results$prcp, method = 'pearson')$estimate, 
  Elevation=cor.test(trend_results$slope, trend_results$Elevation, method = 'pearson')$estimate, 
  Slope=cor.test(trend_results$slope, trend_results$Slope, method = 'pearson')$estimate,
  Aspect=cor.test(trend_results$slope, trend_results$Aspect, method = 'pearson')$estimate,
  Latitude=cor.test(trend_results$slope, trend_results$Latitude, method = 'pearson')$estimate 
)
trend_row[,2:ncol(trend_row)] <- round(trend_row[,2:ncol(trend_row)], 3)

fwrite(rbind(ratio_row, trend_row), paste0(home, "/table_ndvi_cors.csv"))







######################################################################################
## mapping climate correlations and slopes 
results <- fread(paste0(home, "/Data/Climate_NDVI_linmod_r2.csv"))
results$Group <- paste0(results$Variable, " ", results$Season)
results$Group <- as.factor(results$Group)
results$NHDPlusID <- as.numeric(results$NHDPlusID)
results_sig <- results[results$Pvalue <= 0.05,]
results_sig <- results_sig[!results_sig$Slope < quantile(results_sig$Slope, 0.01), ] # get rid of outliers for plotting 
results_sig <- results_sig[!results_sig$Slope > quantile(results_sig$Slope, 0.99), ]

catch <- st_read(paste0(home, "/Data/Catchments/Headwater/headwater_catchments_perm_forest_32617.shp"))
catch <- catch[,c("NHDPlusID", "geometry")]
roi <- st_read("/Volumes/GoogleDrive/My Drive/Chapter2_mechanisms_forest_water_cycling/Data/ROI/blue_ridge_plus_reference.shp")
roi <- st_transform(roi, crs="EPSG:32617")

plot_clim_trends <- function(grp, var, title){
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
  
  gp <- ggplot() + 
    geom_sf(data=roi, fill='white', color='black')+ 
    geom_point(data=grp_cors, aes(x=X, y=Y, color=get(var)), alpha=0.4, size=0.01) + 
    scale_color_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0, name=var) +
    theme_classic() + 
    xlab("") + 
    ylab("") + 
    ggtitle(title) +
    theme(legend.title = element_blank()) + 
    theme(legend.position = "none") +
    #theme(legend.position = c(0.9, 0.3)) + 
    inset_element(trend_violin_plot, 0.5, 0.07, 1, 0.40, align_to = 'full') #left bottom right top
  return(gp)
}


tmin_slope <- plot_clim_trends("tmin WY", "Slope", "Tmin")
tmax_slope <- plot_clim_trends("tmax WY", "Slope", "Tmax")
vp_slope <- plot_clim_trends("vp WY", "Slope", "VP")
prcp_slope <- plot_clim_trends("prcp WY", "Slope", "Precipitation")

WY_results_sig <- results_sig[grep("WY", results_sig$Group),]
WY_results_sig <- WY_results_sig[!WY_results_sig$Variable == "drydays",]
WY_results_sig <- merge(catch, WY_results_sig, inner=TRUE, by="NHDPlusID")
WY_results_sig <- cbind(WY_results_sig, st_coordinates(st_centroid(WY_results_sig)))
legend <- ggplot() + 
  geom_sf(data=roi, fill='white', color='black')+ 
  geom_point(data=WY_results_sig, aes(x=X, y=Y, color=Slope), alpha=0.4, size=0.01) + 
  scale_color_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0, name="tmin WY") +
  theme_classic() + 
  xlab("") + 
  ylab("") + 
  ggtitle("") + 
  theme(legend.title = element_blank()) + 
  theme(legend.position = "right", 
        legend.text=element_text(size=12)) +
  theme(axis.text=element_text(size=12, color = 'black'),
        axis.title=element_text(size=12,face="bold")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1)) 
legend <- ggpubr::get_legend(legend)

tiff(paste0(home, "/Figures/slope_ndvi_clim_models.tiff"), units = 'in', res=200, width = 10, height = 8)
ggarrange(tmin_slope, tmax_slope, vp_slope, prcp_slope,nrow = 2, ncol = 2, 
          common.legend = TRUE,
          legend.grob=legend,
          legend = "right")
dev.off()




plot_clim_trends_rs <- function(grp, var, title){
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
  
  gp <- ggplot() + 
    geom_sf(data=roi, fill='white', color='black')+ 
    geom_point(data=grp_cors, aes(x=X, y=Y, color=get(var)), alpha=0.4, size=0.01) + 
    #scale_color_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0, name=var) +
    scale_color_viridis_c() + 
    theme_classic() + 
    xlab("") + 
    ylab("") + 
    ggtitle(title) +
    theme(legend.title = element_blank()) + 
    theme(legend.position = "none") +
    #theme(legend.position = c(0.9, 0.3)) + 
    inset_element(trend_violin_plot, 0.5, 0.07, 1, 0.40, align_to = 'full') #left bottom right top
  return(gp)
}


tmin_rs <- plot_clim_trends_rs("tmin WY", "Rsquared", "Tmin")
tmax_rs <- plot_clim_trends_rs("tmax WY", "Rsquared", "Tmax")
vp_rs <- plot_clim_trends_rs("vp WY", "Rsquared", "VP")
prcp_rs <- plot_clim_trends_rs("prcp WY", "Rsquared", "Precipitation")

WY_results_sig <- results_sig[grep("WY", results_sig$Group),]
WY_results_sig <- WY_results_sig[!WY_results_sig$Variable == "drydays",]
WY_results_sig <- merge(catch, WY_results_sig, inner=TRUE, by="NHDPlusID")
WY_results_sig <- cbind(WY_results_sig, st_coordinates(st_centroid(WY_results_sig)))
legend <- ggplot() + 
  geom_sf(data=roi, fill='white', color='black')+ 
  geom_point(data=WY_results_sig, aes(x=X, y=Y, color=Rsquared), alpha=0.4, size=0.01) + 
  #scale_color_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0, name="tmin WY") +
  scale_color_viridis_c() + 
  theme_classic() + 
  xlab("") + 
  ylab("") + 
  ggtitle("") + 
  theme(legend.title = element_blank()) + 
  theme(legend.position = "right", 
        legend.text=element_text(size=12)) +
  theme(axis.text=element_text(size=12, color = 'black'),
        axis.title=element_text(size=12,face="bold")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1)) 
legend <- ggpubr::get_legend(legend)

tiff(paste0(home, "/Figures/rsquared_ndvi_clim_models.tiff"), units = 'in', res=200, width = 10, height = 8)
ggarrange(tmin_rs, tmax_rs, vp_rs, prcp_rs,nrow = 2, ncol = 2, 
          common.legend = TRUE,
          legend.grob=legend,
          legend = "right")
dev.off()


# make a table of the quantiles of Rsquared and slope and include fraction of catchments with sig relationship

# Groups at the top: Slope, RSquared, Fraction 
# columns within groups: 1st quartile, median, mean, third quartile, but for fraction it's just the fraction 
# rows are the tmin, tmax, ect
wy_results_sig_dt <- as.data.table(st_drop_geometry(WY_results_sig))
wy_results_sig_dt_summary <- wy_results_sig_dt[, .(fraction_sig = .N/30044, 
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
wy_results_sig_dt_summary[,2:ncol(wy_results_sig_dt_summary)] <- round(wy_results_sig_dt_summary[,2:ncol(wy_results_sig_dt_summary)], 3)
fwrite(wy_results_sig_dt_summary, paste0(home, "/clim_ndvi_model_summary.csv"))



