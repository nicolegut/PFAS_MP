########################################
#### Master's Project Duke MEM 2026
#### UCMR PFAS Interpolations
#### Nicole Gutkowski

##Importing Packages
import pandas as pd
import arcpy, os, time


# time, os, arcpy packages? - used in Hazen

# UCMR Data
path = "C:/Duke/Year 2/MP/Data/Initial Data/ucmr5_occurrence_data"
UCMR_All = "UCMR5_All"
ZIPCodes = "UCMR5_ZIPCodes"
UCMR_data = pd.read_csv(f"{path}/{UCMR_All}.txt", sep="\t", encoding="latin1")
ZIP_data = pd.read_csv(f"{path}/{ZIPCodes}.txt", sep="\t", encoding="latin1")
CWS_data = C:\Duke\Year 2\MP\PFAS_MP.gdb\CWS_Points


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
# look at the duplicate values for pwsid and zip codes
Combine_UCMR = Filtered_UCMR.merge(ZIP_data, on="PWSID", how="left")

UCMR_narm = Combine_UCMR.dropna(subset=["AnalyticalResultValue"])
# checked for duplicates, no duplicate values


##### Breaking the data into compound/ source type of interest
## think this could be turned into some sort of function --> filter(contaminant, source)
## will create 4 dataframes

# subset for contaminant/ water type and check for duplicate combinations

UCMR_PFOA_GW = UCMR_narm[
    (UCMR_narm["Contaminant"] == "PFOA") & (UCMR_narm["FacilityWaterType"] == "GW")
]

UCMR_PFOA_SW = UCMR_narm[
    (UCMR_narm["Contaminant"] == "PFOA") & (UCMR_narm["FacilityWaterType"] == "SW")
]
UCMR_PFOS_GW = UCMR_narm[
    (UCMR_narm["Contaminant"] == "PFOS") & (UCMR_narm["FacilityWaterType"] == "GW")
]
UCMR_PFOS_SW = UCMR_narm[
    (UCMR_narm["Contaminant"] == "PFOS") & (UCMR_narm["FacilityWaterType"] == "SW")
]


# making tuples for combinations of pwsid/zip for each subset
# checking if there are duplicate tuples
# tuple test
PFOA_GW_tuples = [
    tuple(row)
    for row in UCMR_PFOA_GW[["PWSID", "AnalyticalResultValue", "ZIPCODE"]].values
]

"""
add into tuples the concentration data, 
to see if the duplicates have differing concentrations

is there value in keeping the pwsid/ zip uniques? i deleted but could put back
"""

# splitting duplicates
PFOA_GW_dict = {}


# defining a function to assess duplicate PWSID and ZIPCodes for our subset datasets
def count_duplicates(dict, pairs):
    for i in pairs:
        if i in dict:
            dict[i] += 1
        else:
            dict[i] = 1


count_duplicates(PFOA_GW_dict, PFOA_GW_tuples)
"""
okay i have the distinct combinations of zip codes/pwsid/concentration data
need to figure out how many zip codes/ which ones have differing values
- could also be dates? do the same zips get sampled multiple times 

-okay it is dates/ sample event codes - for PWSIDs that sample more than once in the compliance
- determine how many of these there are?

"""


##### Subset the data into test/ control groups
## probably still in pandas
## CONSIDER: for subset - look into what is a good proportion to subset/ even spatial sampling
##will create 8 dataframes


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
