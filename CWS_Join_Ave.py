### CWS Join Script
###Goal: to make four datasets of points


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


basepath = "C:/Duke/Year 2/MP/"
fgdb = os.path.join(basepath, "PFAS_MP.gdb")
CWS_pts = os.path.join(fgdb, "CWS_Points")


## Zip code data need to make sure they have leading zeros!!
## check both UCMR and Hazen data


# Hazen Data
# zipcode data
# points of zipcode centroid points (contiguous US)
# shapefile of all contiguous US zipcodes
# line around country (dissolved from shp - contiguous US)
# mask of counties - raster including AK + HI


#######################################################################
############## Thought Process/ Plan ##################################

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

# think about this --> also throws a warning about copies/ slices
Filtered_UCMR["AnalyticalResultValue"] = Filtered_UCMR["AnalyticalResultValue"].fillna(
    0.004 * (1 / 3)
)


def df_to_gdb(df, gdb_path, table_name):
    # Create a table path
    table_path = f"{gdb_path}/{table_name}"

    # Create a table using arcpy
    arcpy.management.CreateTable(gdb_path, table_name)

    # Dict of Field Types
    fieldtypes = {
        "datetime64[ns]": "DATE",
        "int64": "LONG",
        "object": "TEXT",
        "float64": "DOUBLE",
    }

    # Add fields based on DataFrame columns
    for col in df.columns:
        arcpy.management.AddField(table_path, col, fieldtypes[df[col].dtype.name])

    # Insert data from DataFrame to the GDB table
    cursor = arcpy.da.InsertCursor(table_path, df.columns.tolist())
    for index, row in df.iterrows():
        cursor.insertRow(tuple(row))
    del cursor

    # return path
    return table_path


##Analysis

# edit these: ryan's code for setting workspaces/ environments

print("Setting Workspace")
arcpy.env.workspace = fgdb
arcpy.env.overwriteOutput = True
# arcpy.env.parallelProcessingFactor = "100%"


##create the dataframes for each and take the max value per PWSID
# only take the max after separating the dataframes

UCMR_PFOA_GW = Filtered_UCMR[
    (Filtered_UCMR["Contaminant"] == "PFOA")
    & (Filtered_UCMR["FacilityWaterType"] == "GW")
]

UCMR_PFOA_SW = Filtered_UCMR[
    (Filtered_UCMR["Contaminant"] == "PFOA")
    & (Filtered_UCMR["FacilityWaterType"] == "SW")
]

UCMR_PFOS_GW = Filtered_UCMR[
    (Filtered_UCMR["Contaminant"] == "PFOS")
    & (Filtered_UCMR["FacilityWaterType"] == "GW")
]
UCMR_PFOS_SW = Filtered_UCMR[
    (Filtered_UCMR["Contaminant"] == "PFOS")
    & (Filtered_UCMR["FacilityWaterType"] == "SW")
]


def PWSID_summary_stats_calc(df):
    # grouping by PWSID
    group_cols = "PWSID"
    # taking summary stats of analytical result value
    metric_cols = "AnalyticalResultValue"
    # columns that we would like to maintain in the dataframes
    keep_cols = [
        "PWSID",
        "PWSName",
        "Size",
        "FacilityWaterType",
        "Contaminant",
        "MRL",
        "Units",
        "Region",
        "State",
    ]

    # creating a dataframe with only the columns that we could like to keep
    # dropping duplicate rows in the columns that we wanted to maintain
    # (one obs per PWSID + adding on the summary stats after)
    df1 = df[keep_cols].drop_duplicates(subset=group_cols, keep="first").copy()

    # summary stats!
    # mean
    means = df.groupby(group_cols)[metric_cols].mean()
    # min
    mins = df.groupby(group_cols)[metric_cols].min()
    # max
    maxes = df.groupby(group_cols)[metric_cols].max()
    # range
    ranges = maxes - mins
    # stdev
    stdev = df.groupby(group_cols)[metric_cols].std()

    # merging/ renaming the summary stats in a table
    merge_stats = [
        means.rename("MeanValue"),
        mins.rename("MinValue"),
        maxes.rename("MaxValue"),
        ranges.rename("Range"),
        stdev.rename("StdDev"),
    ]
    merge_stats = pd.concat(merge_stats, axis=1)

    # merging the summary stats table by PWSID to the cleaned df
    final_df = df1.merge(
        right=merge_stats, right_index=True, left_on=group_cols, how="right"
    ).copy()
    return final_df


PFOA_GW_stats = PWSID_summary_stats_calc(UCMR_PFOA_GW)
PFOA_SW_stats = PWSID_summary_stats_calc(UCMR_PFOA_SW)
PFOS_GW_stats = PWSID_summary_stats_calc(UCMR_PFOS_GW)
PFOS_SW_stats = PWSID_summary_stats_calc(UCMR_PFOS_SW)


pfoa_sw_tbl = df_to_gdb(PFOA_SW_stats, fgdb, "PFOA_SW")
pfoa_gw_tbl = df_to_gdb(PFOA_GW_stats, fgdb, "PFOA_GW")
pfos_sw_tbl = df_to_gdb(PFOS_SW_stats, fgdb, "PFOS_SW")
pfos_gw_tbl = df_to_gdb(PFOS_GW_stats, fgdb, "PFOS_GW")

## Join to Zips

table_list = [
    pfoa_gw_tbl,
    pfos_sw_tbl,
    pfos_gw_tbl,
    pfoa_sw_tbl,
]  # This list can be altered to delineate IDWs for select contaminants
FC_List = []


for table in table_list:
    table_name = str(table).replace("_tbl", "_output").replace("//", "/")
    print("table_name is " + table_name)
    outpath = table_name + "_Oct2025"
    print("outpath is " + outpath)
    if os.path.exists(outpath):
        os.remove(outpath)
    joinedtable = arcpy.management.AddJoin(
        CWS_pts, "PWSID", table, "PWSID", "KEEP_COMMON"
    )
    print(joinedtable)
    arcpy.management.CopyFeatures(joinedtable, outpath)
    FC_List.append(outpath)

## Interpolate with IDW
# Set environment settings
print(FC_List)
arcpy.env.workspace = fgdb
