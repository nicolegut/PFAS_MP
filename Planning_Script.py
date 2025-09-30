########################################
#### Master's Project Duke MEM 2026
#### UCMR PFAS Interpolations
#### Nicole Gutkowski

##Importing Packages
import pandas as pd

# time, os, arcpy packages? - used in Hazen

# UCMR Data
path = "C:/Duke/Year 2/MP/Data/Initial Data/ucmr5_occurrence_data"
UCMR_All = "UCMR5_All"
ZIPCodes = "UCMR5_ZIPCodes"
UCMR_data = pd.read_csv(f"{path}/{UCMR_All}.txt", sep="\t", encoding="latin1")
ZIP_data = pd.read_csv(f"{path}/{ZIPCodes}.txt", sep="\t", encoding="latin1")


## Zip code data need to make sure they have leading zeros!!
## check both UCMR and Hazen data


# Hazen Data
# zipcode data
# points of zipcode centroid points (contiguous US)
# shapefile of all contiguous US zipcodes
# line around country (dissolved from shp - contiguous US)
# mask of counties - raster including AK + HI


############################################################################
############### Thought Process/ Plan #####################################

#### Data Wrangling: ##################################################


##### Upload UCMR and PWSID+zip data

##make sure to put in leading zeros
ZIP_data[["ZIPCODE"]] = ZIP_data[["ZIPCODE"]].astype(str)
ZIP_data[["ZIPCODE"]] = ZIP_data[["ZIPCODE"]].apply(lambda x: x.str.zfill(5))

##filter out unwanted contaminants and water types
Filtered_UCMR = UCMR_data[
    (UCMR_data["Contaminant"].isin(["PFOA", "PFOS"]))
    & (UCMR_data["FacilityWaterType"].isin(["GW", "SW"]))
]

## Merge the datasets together (join by zip code? hopefully)
##one facility serves muliple zip codes --> need to figure out how to best do this
# Combine_UCMR = pd.merge(ZIP_data, UCMR_data, on="PWSID", how="XXXXX")

# Merge with ZIP codes
Combine_UCMR = Filtered_UCMR.merge(ZIP_data, on="PWSID", how="left")

# thoughts on the merge


##### Breaking the data into compound/ source type of interest
## think this might be a some sort of function --> filter(contaminant, source)
## will create 4 dataframes

##### Subset the data into test/ control groups
## probably still in pandas
## CONSIDER: for subset - look into what is a good proportion to subset/ even spatial sampling
##will create 8 dataframes
UCMR_PFOA_GW = UCMR_data[
    (UCMR_data["Contaminant"] == "PFOA") & (UCMR_data["FacilityWaterType"] == "GW")
]
# UCMR_PFOA_SW =
# UCMR_PFOS_GW =
# UCMR_PFOS_SW =

#### Output of Data Wrangling:
"""
4 test df: PFOA+SW/ PFOA+GW/ PFOS+SW/ PFOS+GW
4 check df - subset of the above
"""

#### Interpolation Analysis: #########################################

####converting the dataframes into gdb tables
## in ryan's code --> df_to_gdb
##


# can each interpolation be run as a function?
## input would be the dataset
## for NA's --> in ryans script he filled with 0.004 * (1/3) = 0.00133
## considerations:
#### environment settings
#### groundwater/ surface water differences --> which interpolation methods can actually input the 'barriers'


#### Output of Interpolation Analysis
"""
For each interpolation method
- 4 interpolations
- gdb points used to create interpolation
- gdb points used to check interpolation
    
"""

###### Stat Analysis ###########################
## using the test points
## create new field --> pull value from raster
## will have fields with actual and predicted value
