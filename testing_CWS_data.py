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

"""
first try the 
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

##putting back in the filtered UCMR data
Filtered_UCMR = UCMR_data[
    (UCMR_data["Contaminant"].isin(["PFOA", "PFOS"]))
    & (UCMR_data["FacilityWaterType"].isin(["GW", "SW"]))
]

# finding unique PWSID for the UCMR and CWS data
UCMR_PWSID = set(Filtered_UCMR["PWSID"])
CWS_PWSID = set(CWS_data["PWSID"])

# finding how many of the UCMR pwsids are represented in the CWS pwsids
# left merge, so non-matches in the CWS dataset would be blank
PWSID_check = UCMR_PWSID & CWS_PWSID
unaccounted_PWSID = len(PWSID_check) - len(UCMR_PWSID)
print(unaccounted_PWSID)

##compare against how many pwsid are served by multiple zip codes
##sum the 'data losses/ how many PWSIDs we would lose by choosing 1 per zip
