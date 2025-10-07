##testing the use of the CWS data rather than zipcodes
# how many would we be losing in comparison to using the zip data

"""
# of PWSID data that would be unaccounted for using the zip codes
versus
PWSID unaccounted for using the CWS data
"""

"""
thinking I do the same tuple thing, find the unique combinations 
and then create a new dict that should how many pwsids are assoc with the 
same zip codes --> the difference between the unique tuples and the # of values in the
new zip dictionary would be how many pwsid unique datapoints would be unaccounted for
"""


import pandas as pd
import arcpy, os, time

# time, os, arcpy packages? - used in Hazen

# UCMR Data
path = "C:/Duke/Year 2/MP/Data/Initial Data/ucmr5_occurrence_data"
UCMR_All = "UCMR5_All"
ZIPCodes = "UCMR5_ZIPCodes"
UCMR_data = pd.read_csv(f"{path}/{UCMR_All}.txt", sep="\t", encoding="latin1")
ZIP_data = pd.read_csv(f"{path}/{ZIPCodes}.txt", sep="\t", encoding="latin1")
CWS = "C:/Duke/Year 2/MP/PFAS_MP.gdb/CWS_Points"

##filling in the zip data with leading 0's
ZIP_data[["ZIPCODE"]] = ZIP_data[["ZIPCODE"]].astype(str)
ZIP_data[["ZIPCODE"]] = ZIP_data[["ZIPCODE"]].apply(lambda x: x.str.zfill(5))


##not my function --> from d-wasserman/FeatureTabletoDataframe.py on github gist!
##saved in pythonref in edge favorites folder
def arcgis_table_to_df(in_fc, input_fields=None, query=""):
    """Function will convert an arcgis table into a pandas dataframe with an object ID index, and the selected
    input fields using an arcpy.da.SearchCursor.
    :param - in_fc - input feature class or table to convert
    :param - input_fields - fields to input to a da search cursor for retrieval
    :param - query - sql query to grab appropriate values
    :returns - pandas.DataFrame"""
    OIDFieldName = arcpy.Describe(in_fc).OIDFieldName
    if input_fields:
        final_fields = [OIDFieldName] + input_fields
    else:
        final_fields = [field.name for field in arcpy.ListFields(in_fc)]
    data = [
        row for row in arcpy.da.SearchCursor(in_fc, final_fields, where_clause=query)
    ]
    fc_dataframe = pd.DataFrame(data, columns=final_fields)
    fc_dataframe = fc_dataframe.set_index(OIDFieldName, drop=True)
    return fc_dataframe


# getting the CWS point values into a pandas df
CWS_data = arcgis_table_to_df(CWS)

"""
These steps:
- filter the data into contaminants (PFOA and PFOS) and water types (GW and SW)
- 
"""

##putting back in the filtered UCMR data
Filtered_UCMR = UCMR_data[
    (UCMR_data["Contaminant"].isin(["PFOA", "PFOS"]))
    & (UCMR_data["FacilityWaterType"].isin(["GW", "SW"]))
]

Combine_UCMR = Filtered_UCMR.merge(ZIP_data, on="PWSID", how="left")


Zip_PWSID_tuples = [tuple(row) for row in Combine_UCMR[["PWSID", "ZIPCODE"]].values]

mult_PWSID_to_ZIP_dict = {}


def count_PWSID_serving_ZIPs(dict, pairs):
    for pwsid, zipcode in pairs:
        if zipcode not in dict:
            dict[zipcode] = {pwsid}
        else:
            dict[zipcode].add(pwsid)


count_PWSID_serving_ZIPs(mult_PWSID_to_ZIP_dict, Zip_PWSID_tuples)

## ^^^^ investigate this further? idk


# finding unique PWSID for the UCMR and CWS data
UCMR_PWSID = set(Filtered_UCMR["PWSID"])
CWS_PWSID = set(CWS_data["PWSID"])

# finding how many of the UCMR pwsids are represented in the CWS pwsids
# left merge, so non-matches in the CWS dataset would be blank
PWSID_check = UCMR_PWSID & CWS_PWSID
unaccounted_PWSID = len(UCMR_PWSID) - len(PWSID_check)
print(unaccounted_PWSID)
##448 unnacounted for PWSIDs using the CWS data as a spatial match

##seeing if all of the PWSIDs are represented by the zip code data anyways
ZIP_PWSID = set(ZIP_data["PWSID"])
PWSID_check2 = UCMR_PWSID & ZIP_PWSID
unaccounted_PWSID2 = len(UCMR_PWSID) - len(PWSID_check2)
## 146 unrepresented PWSIDs using zip codes as a spatial match

##checking how many


##compare against how many pwsid are served by multiple zip codes
##sum the 'data losses/ how many PWSIDs we would lose by choosing 1 per zip

UCMR_PFOA_GW = Combine_UCMR[
    (Combine_UCMR["Contaminant"] == "PFOA")
    & (Combine_UCMR["FacilityWaterType"] == "GW")
]
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


# making tuples for combinations of pwsid/zip for each subset
# checking if there are duplicate tuples
# tuple test
PFOA_GW_tuples = [
    tuple(row)
    for row in UCMR_PFOA_GW[["PWSID", "AnalyticalResultValue", "ZIPCODE"]].values
]
PFOA_SW_tuples = [
    tuple(row)
    for row in UCMR_PFOA_SW[["PWSID", "AnalyticalResultValue", "ZIPCODE"]].values
]
PFOS_GW_tuples = [
    tuple(row)
    for row in UCMR_PFOS_GW[["PWSID", "AnalyticalResultValue", "ZIPCODE"]].values
]
PFOS_SW_tuples = [
    tuple(row)
    for row in UCMR_PFOS_SW[["PWSID", "AnalyticalResultValue", "ZIPCODE"]].values
]

PFOA_GW_dict = {}
PFOA_SW_dict = {}
PFOS_GW_dict = {}
PFOS_SW_dict = {}


# defining a function to assess duplicate PWSID and ZIPCodes for our subset datasets
def count_duplicates(dict, pairs):
    for i in pairs:
        if i in dict:
            dict[i] += 1
        else:
            dict[i] = 1


count_duplicates(PFOA_GW_dict, PFOA_GW_tuples)
count_duplicates(PFOA_SW_dict, PFOA_SW_tuples)
count_duplicates(PFOS_GW_dict, PFOS_GW_tuples)
count_duplicates(PFOS_SW_dict, PFOS_SW_tuples)

# counting how many
PFOA_GW_pwsidict = {}
PFOA_SW_pwsiddict = {}
PFOS_GW_pwsiddict = {}
PFOS_SW_pwsiddict = {}


##checking mult pwsids in original zip data
