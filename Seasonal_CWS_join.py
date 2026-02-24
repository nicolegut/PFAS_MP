### CWS Join Script - Seasonal
###Goal: to make four datasets of points
##use arcpy env to run this one!

##Importing Packages
import pandas as pd
import arcpy, os, time


# UCMR Data
path = "C:/Duke/Year 2/MP/Data/Initial Data/ucmr5_occurrence_data"
UCMR_All = "UCMR5_All"
UCMR_data = pd.read_csv(f"{path}/{UCMR_All}.txt", sep="\t", encoding="latin1")
CWS_data = "C:/Duke/Year 2/MP/PFAS_MP.gdb/CWS_Points"


basepath = "C:/Duke/Year 2/MP/"
fgdb = os.path.join(basepath, "PFAS_MP.gdb")
# CWS points were created by the 'points from polygon centroid' in arcpro/ not in a script
# centroid points were forced to be inside the boundary of their polygons
CWS_pts = os.path.join(fgdb, "CWS_Points")
seas_gdb = os.path.join(basepath, "Seasonal_Tests.gdb")
us_mask = os.path.join(fgdb, "ContigUS_Mask")


# Hazen Data
# zipcode data
# points of zipcode centroid points (contiguous US)
# shapefile of all contiguous US zipcodes
# line around country (dissolved from shp - contiguous US)
# mask of counties - raster including AK + HI


#######################################################################
############## Thought Process/ Plan ##################################

#### Data Wrangling: ##################################################


##### Filter UCMR data

##filter out unwanted contaminants and water types
Filtered_UCMR = UCMR_data[
    (UCMR_data["Contaminant"].isin(["PFOA", "PFOS"]))
    & (UCMR_data["FacilityWaterType"].isin(["GW", "SW"]))
]

# think about this --> also throws a warning about copies/ slices
Filtered_UCMR["AnalyticalResultValue"] = Filtered_UCMR["AnalyticalResultValue"].fillna(
    0.004 * (1 / 3)
)

Filtered_UCMR["CollectionDate"] = pd.to_datetime(Filtered_UCMR["CollectionDate"])
Filtered_UCMR["Month"] = Filtered_UCMR["CollectionDate"].dt.month


# function for converting df to gdb dataframe
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


summer_months = [6, 7, 8]
winter_months = [12, 1, 2]

##Analysis

print("Setting Workspace")
arcpy.env.workspace = fgdb
arcpy.env.overwriteOutput = True
# arcpy.env.parallelProcessingFactor = "100%"


##create the dataframes for each contaminant/ water pair
# PFOA GW
PFOA_GW_sum = Filtered_UCMR[
    (Filtered_UCMR["Contaminant"] == "PFOA")
    & (Filtered_UCMR["FacilityWaterType"] == "GW")
    & (Filtered_UCMR["Month"].isin(summer_months))
]

PFOA_GW_win = Filtered_UCMR[
    (Filtered_UCMR["Contaminant"] == "PFOA")
    & (Filtered_UCMR["FacilityWaterType"] == "GW")
    & (Filtered_UCMR["Month"].isin(winter_months))
]

# PFOA SW
PFOA_SW_sum = Filtered_UCMR[
    (Filtered_UCMR["Contaminant"] == "PFOA")
    & (Filtered_UCMR["FacilityWaterType"] == "SW")
    & (Filtered_UCMR["Month"].isin(summer_months))
]

PFOA_SW_win = Filtered_UCMR[
    (Filtered_UCMR["Contaminant"] == "PFOA")
    & (Filtered_UCMR["FacilityWaterType"] == "SW")
    & (Filtered_UCMR["Month"].isin(winter_months))
]

# PFOS_GW
PFOS_GW_sum = Filtered_UCMR[
    (Filtered_UCMR["Contaminant"] == "PFOS")
    & (Filtered_UCMR["FacilityWaterType"] == "GW")
    & (Filtered_UCMR["Month"].isin(summer_months))
]

PFOS_GW_win = Filtered_UCMR[
    (Filtered_UCMR["Contaminant"] == "PFOS")
    & (Filtered_UCMR["FacilityWaterType"] == "GW")
    & (Filtered_UCMR["Month"].isin(winter_months))
]

# PFOS_SW
PFOS_SW_sum = Filtered_UCMR[
    (Filtered_UCMR["Contaminant"] == "PFOS")
    & (Filtered_UCMR["FacilityWaterType"] == "SW")
    & (Filtered_UCMR["Month"].isin(summer_months))
]

PFOS_SW_win = Filtered_UCMR[
    (Filtered_UCMR["Contaminant"] == "PFOS")
    & (Filtered_UCMR["FacilityWaterType"] == "SW")
    & (Filtered_UCMR["Month"].isin(winter_months))
]


# function for adding summary stats to the filtered dfs
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

    df["CollectionDate"] = pd.to_datetime(df["CollectionDate"], format="%m/%d/%Y")

    df1 = df[keep_cols].drop_duplicates(subset=group_cols, keep="first").copy()

    # summary stats!
    # mean
    means = df.groupby(group_cols)[metric_cols].mean()
    mins = df.groupby(group_cols)[metric_cols].min()
    # max
    maxes = df.groupby(group_cols)[metric_cols].max()
    # range
    ranges = maxes - mins
    # stdev
    stdev = df.groupby(group_cols)[metric_cols].std()
    # count
    counts = df.groupby(group_cols)[metric_cols].count()
    # date ranges
    min_date = df.groupby("PWSID")["CollectionDate"].min()
    max_date = df.groupby("PWSID")["CollectionDate"].max()

    date_count = df.groupby("PWSID")["CollectionDate"].nunique()

    # merging/ renaming the summary stats in a table
    merge_stats = [
        means.rename("MeanValue"),
        mins.rename("MinValue"),
        maxes.rename("MaxValue"),
        ranges.rename("Range"),
        stdev.rename("StdDev"),
        counts.rename("CountofPoints"),
        date_count.rename("CountofDates"),
        min_date.rename("FirstCollectionDate"),
        max_date.rename("LastCollectionDate"),
    ]

    merge_stats = pd.concat(merge_stats, axis=1)

    # merging the summary stats table by PWSID to the cleaned df
    final_df = df1.merge(
        right=merge_stats, right_index=True, left_on=group_cols, how="right"
    ).copy()
    return final_df


# adding summary stats to the filtered datasets
PFOA_GW_sum_stats = PWSID_summary_stats_calc(PFOA_GW_sum)
PFOA_SW_sum_stats = PWSID_summary_stats_calc(PFOA_SW_sum)
PFOS_GW_sum_stats = PWSID_summary_stats_calc(PFOS_GW_sum)
PFOS_SW_sum_stats = PWSID_summary_stats_calc(PFOS_SW_sum)

PFOA_GW_win_stats = PWSID_summary_stats_calc(PFOA_GW_win)
PFOA_SW_win_stats = PWSID_summary_stats_calc(PFOA_SW_win)
PFOS_GW_win_stats = PWSID_summary_stats_calc(PFOS_GW_win)
PFOS_SW_win_stats = PWSID_summary_stats_calc(PFOS_SW_win)


# exporting to gdb for use in arc
pfoa_sw_sum = df_to_gdb(PFOA_SW_sum_stats, seas_gdb, "PFOA_SW_sum")
pfoa_gw_sum = df_to_gdb(PFOA_GW_sum_stats, seas_gdb, "PFOA_GW_sum")
pfos_sw_sum = df_to_gdb(PFOS_SW_sum_stats, seas_gdb, "PFOS_SW_sum")
pfos_gw_sum = df_to_gdb(PFOS_GW_sum_stats, seas_gdb, "PFOS_GW_sum")

pfoa_sw_win = df_to_gdb(PFOA_SW_win_stats, seas_gdb, "PFOA_SW_win")
pfoa_gw_win = df_to_gdb(PFOA_GW_win_stats, seas_gdb, "PFOA_GW_win")
pfos_sw_win = df_to_gdb(PFOS_SW_win_stats, seas_gdb, "PFOS_SW_win")
pfos_gw_win = df_to_gdb(PFOS_GW_win_stats, seas_gdb, "PFOS_GW_win")


## Join to CWS points
table_list = [
    pfoa_gw_sum,
    pfos_sw_sum,
    pfos_gw_sum,
    pfoa_sw_sum,
    pfoa_gw_win,
    pfos_sw_win,
    pfos_gw_win,
    pfoa_sw_win,
]
FC_List = []

# joining the CWS points to the gdb tables - need to change the outpath name for future runs
for table in table_list:
    base_name = os.path.basename(table)
    outpath = os.path.join(seas_gdb, f"{base_name}_Feb25")
    print(f"table_name is {base_name}"),
    print(f"outpath {outpath}"),
    if arcpy.Exists(outpath):
        arcpy.management.Delete(outpath)
    joinedtable = arcpy.management.AddJoin(
        CWS_pts, "PWSID", table, "PWSID", "KEEP_COMMON"
    )
    print(joinedtable)
    arcpy.management.CopyFeatures(joinedtable, outpath)
    FC_List.append(outpath)


PFOAGW_win = os.path.join(seas_gdb, "PFOA_GW_win_Feb25")
PFOASW_win = os.path.join(seas_gdb, "PFOA_SW_win_Feb25")
PFOSGW_win = os.path.join(seas_gdb, "PFOS_GW_win_Feb25")
PFOSSW_win = os.path.join(seas_gdb, "PFOS_SW_win_Feb25")

PFOAGW_sum = os.path.join(seas_gdb, "PFOA_GW_sum_Feb25")
PFOASW_sum = os.path.join(seas_gdb, "PFOA_SW_sum_Feb25")
PFOSGW_sum = os.path.join(seas_gdb, "PFOS_GW_sum_Feb25")
PFOSSW_sum = os.path.join(seas_gdb, "PFOS_SW_sum_Feb25")

point_list = [
    PFOAGW_win,
    PFOSGW_win,
    PFOASW_win,
    PFOSSW_win,
    PFOAGW_sum,
    PFOSGW_sum,
    PFOASW_sum,
    PFOSSW_sum,
]

naming_list = [
    "PFOA_GW_win",
    "PFOS_GW_win",
    "PFOA_SW_win",
    "PFOS_SW_win",
    "PFOA_GW_sum",
    "PFOS_GW_sum",
    "PFOA_SW_sum",
    "PFOS_SW_sum",
]

for point, name in zip(point_list, naming_list):
    out_train = os.path.join(seas_gdb, f"{name}_training")
    out_test = os.path.join(seas_gdb, f"{name}_testing")
    contig_pts = os.path.join(seas_gdb, f"{name}_contig")

    arcpy.analysis.Clip(
        in_features=point,
        clip_features=us_mask,
        out_feature_class=contig_pts,
        cluster_tolerance=None,
    )

    arcpy.SubsetFeatures_ga(
        in_features=f"{contig_pts}",
        out_training_feature_class=out_train,
        out_test_feature_class=out_test,
        size_of_training_dataset="90",
        subset_size_units="PERCENTAGE_OF_INPUT",
    )
