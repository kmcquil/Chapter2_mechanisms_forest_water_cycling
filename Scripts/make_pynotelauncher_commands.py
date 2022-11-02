from glob import glob
import os 
import geopandas as gpd
import pandas as pd
home = "/share/klmarti3/kmcquil/Chapter2_mechanisms_forest_water_cycling"

# make a list of all of the nhdplus ids from the headwater catchments shp that are permanently forested 
hws = gpd.read_file(glob(os.path.join(home, "Data", "Catchments", "Headwater", "headwater_catchments_perm_forest_32617.shp"))[0])
ids = hws['NHDPlusID']

rows = []
for id in ids:
    rows.append('python /share/klmarti3/kmcquil/Chapter2_mechanisms_forest_water_cycling/Scripts/ratio_ndvi_landsat_mpi.py ' + str(id))

df = pd.DataFrame(rows)
df.to_csv(os.path.join(home, "Scripts", "commands_ratio_ndvi_landsat_mpi_permanent_forest.txt"), header = None, index=None)



# make a list of all of the nhdplus ids from the headwater catchments shp that are 100% forested 
c100 = pd.read_csv(os.path.join(home, "Scripts", "commands_ratio_ndvi_landsat_mpi_permanent_forest.txt"), header=None)
print(c100.shape)

# find which of the IDs was finished 
result_filepaths = glob(os.path.join(home, "Data", "NDVI", "catchment_ratio_ndvi_results", "*.csv"))
ids_done = [s[115:129] for s in result_filepaths]
ids_done = '|'.join(ids_done)

# filter them out of the list 
to_do = c100[c100[0].str.contains(ids_done) == False]
print(to_do.shape)

# save the new list 
to_do.to_csv(os.path.join(home, "Scripts", "commands_ratio_ndvi_landsat_mpi_permanent_forest_NOT_FINISHED.txt"), header = None, index=None)



####################################################################################################################################################################
####################################################################################################################################################################
## Make the commands for the ET mpi tasks 
##rows = []
##for id in ids:
##    rows.append('python /share/klmarti3/kmcquil/Chapter2_mechanisms_forest_water_cycling/Scripts/et_catchment_metrics_mpi.py ' + str(id))

##df = pd.DataFrame(rows)
##df.to_csv(os.path.join(home, "Scripts", "commands_et_metrics_mpi_permanent_forest.txt"), header = None, index=None)

# make a list of all of the nhdplus ids from the headwater catchments shp that are 100% forested 
##c100 = pd.read_csv(os.path.join(home, "Scripts", "commands_et_metrics_mpi_permanent_forest.txt"), header=None)
##print(c100.shape)

# find which of the IDs was finished 
##result_filepaths = glob(os.path.join(home, "Data", "Evapotranspiration", "catchment_et_metrics", "*.csv"))
##ids_done = [s[119:133] for s in result_filepaths]
##ids_done = '|'.join(ids_done)

# filter them out of the list 
##to_do = c100[c100[0].str.contains(ids_done) == False]
##print(to_do.shape)

##to_do.to_csv(os.path.join(home, "Scripts", "commands_et_metrics_mpi_permanent_forest_NOT_FINISHED.txt"), header = None, index=None)






####################################################################################################################################################################
####################################################################################################################################################################
## Check which catchments ET didn't get summarized 

##et = pd.read_csv(os.path.join(home, "Data","Evapotranspiration","summary_et.csv"))

# make a list of all of the nhdplus ids from the headwater catchments shp that are 100% forested 
##c100 = pd.read_csv(os.path.join(home, "Scripts", "commands_et_metrics_mpi_permanent_forest.txt"), header=None)
##print(c100.shape)

# find the IDs that weren't messed up 
##et.wsid = et.wsid.apply(str)
##ids_done = '|'.join(et.wsid)

# filter them out of the list 
##to_do = c100[c100[0].str.contains(ids_done) == False]
##print(to_do.shape)

##to_do.to_csv(os.path.join(home, "Scripts", "commands_et_metrics_mpi_permanent_forest_NOT_FINISHED.txt"), header = None, index=None)

