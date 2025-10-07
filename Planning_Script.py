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
CWS_data = "C:/Duke/Year 2/MP/PFAS_MP.gdb/CWS_Points"


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

# Merge with ZIP codes
# look at the duplicate values for pwsid and zip codes
Combine_UCMR = Filtered_UCMR.merge(ZIP_data, on="PWSID", how="left")

##almost 5.5x more rows when joined because of duplicate dates/ zipcodes associated with the PWSIDs


"""
This code assesses multiple PWSIDs serving single zip codes
the count_duplicates function just provides a count of pwsids associated with zips
the count_PWSIDs_serving_zips function provides the names of the PWSIDs
"""

PWSID_Zip_dupl = {}
PWSID_Zip_duplID = {}
PWSID_ZIP_tuples = [tuple(row) for row in ZIP_data[["PWSID", "ZIPCODE"]].values]


# defining a function to assess duplicate PWSID and ZIPCodes for our subset datasets
# This function just provides a count of the pwsids
def count_duplicates(dict, pairs):
    for i, j in pairs:
        if j in dict:
            dict[j] += 1
        else:
            dict[j] = 1


##This function provides the actual PWSID numbers that are associated with the Zip
def count_PWSID_serving_ZIPs(dict, pairs):
    for pwsid, zipcode in pairs:
        if zipcode not in dict:
            dict[zipcode] = {pwsid}
        else:
            dict[zipcode].add(pwsid)


count_duplicates(PWSID_Zip_dupl, PWSID_ZIP_tuples)

count_PWSID_serving_ZIPs(PWSID_Zip_duplID, PWSID_ZIP_tuples)


"""
Using the combined data, potentially taking max value for each 
Ideas: taking the max, min, median or mean 
- median/ mean would lose the PWSID association

group by Zipcodes
- need to maintain other values for max and mins
"""


def max_conc_zip(df, contaminant, watertype):
    filtered_df = df[
        (df["Contaminant"] == contaminant) & (df["FacilityWaterType"] == watertype)
    ]
    grouped_df = filtered_df.groupby("ZIPCODE").max("AnalyticalResultValue")
    return grouped_df


Max_PFOA_GW = max_conc_zip(Combine_UCMR, "PFOA", "GW")


##### Breaking the data into compound/ source type of interest
## think this could be turned into some sort of function --> filter(contaminant, source)
## will create 4 dataframes

# subset for contaminant/ water type and check for duplicate combinations
# Using the current combined data, but may change to

UCMR_PFOA_GW = Combine_UCMR[
    (Combine_UCMR["Contaminant"] == "PFOA")
    & (Combine_UCMR["FacilityWaterType"] == "GW")
]

PFOA_GW_maxtest = UCMR_PFOA_GW.groupby("ZIPCODE").max("AnalyticalResultValue")

UCMR_PFOA_SW = Combine_UCMR[
    (Combine_UCMR["Contaminant"] == "PFOA")
    & (Combine_UCMR["FacilityWaterType"] == "SW")
]
UCMR_PFOS_GW = Combine_UCMR[
    (Combine_UCMR["Contaminant"] == "PFOS")
    & (Combine_UCMR["FacilityWaterType"] == "GW")
]
UCMR_PFOS_SW = Combine_UCMR[
    (Combine_UCMR["Contaminant"] == "PFOS")
    & (Combine_UCMR["FacilityWaterType"] == "SW")
]


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
