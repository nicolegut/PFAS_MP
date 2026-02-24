##Importing Packages
import pandas as pd
import arcpy, os, time
from datetime import datetime
import itertools
import numpy as np
from itertools import cycle


# setting up paths
basepath = "C:/Duke/Year 2/MP/"
fgdb = os.path.join(basepath, "PFAS_MP.gdb")
contig_gdb = os.path.join(fgdb, "Exp_Int_pts")
testing_gdb = os.path.join(fgdb, "Testing_Pts")
# change this one for different interpolations

PFOA_GW_folder = os.path.join(basepath, r"Interpolation_testing\PFOA_GW")


# paths to points for subsetting
train_pts = os.path.join(testing_gdb, "PFOA_GW_training")
test_pts = os.path.join(testing_gdb, "PFOA_GW_testing")
test_bufs5mi = os.path.join(testing_gdb, "PFOA_GW_test_5mibuf")
us_mask = os.path.join(fgdb, "ContigUS_Mask")
us_snap_ras = os.path.join(fgdb, "ContigUS_Raster")

# top 75 rasters from the
PFOA_GW_IDW = os.path.join(PFOA_GW_folder, "IDW_75_Params.csv")
PFOA_GW_UK = os.path.join(PFOA_GW_folder, "UK_75_Params.csv")
PFOA_GW_OK = os.path.join(PFOA_GW_folder, "OK_75_Params.csv")

arcpy.env.overwriteOutput = True
arcpy.env.mask = us_mask
arcpy.env.extent = us_mask
arcpy.env.snapRaster = us_snap_ras
arcpy.env.cellSize = us_snap_ras
# cell size is 5 mi in meters (8046.7m)

print(f"set environments/ loading csv for surface water")

mapping = {"8000": "8", "16000": "16", "20000": "20", "24000": "24", "#": "Arc"}


# read csv for top 75 and edit so that the params match the iterators
IDW_params = pd.read_csv(PFOA_GW_IDW)
IDW_params["out_name"] = IDW_params["raster"].apply(os.path.basename)

OK_params = pd.read_csv(PFOA_GW_OK)
OK_params["out_name"] = OK_params["raster"].apply(os.path.basename)
OK_params["lag_name"] = OK_params["lag_dist"].map(mapping)

UK_params = pd.read_csv(PFOA_GW_UK)
UK_params["out_name"] = UK_params["raster"].apply(os.path.basename)
UK_params["lag_name"] = UK_params["lag_dist"].map(mapping)


# field interpolating
z_field = "MeanValue"

overall_script_start = time.time()


print(
    f"loaded parameter tables/ starting IDW Now [{datetime.now().strftime('%H:%M:%S')}]"
)


########################################################################################
"""
IDW Testing - Calculating Stats
"""
########################################################################################

print("calculating stats for IDW")


## setting up defs for analysis
def get_alias_not_field_name(fc, alias):
    for f in arcpy.ListFields(fc):
        if f.aliasName == alias:
            return f.name
    raise ValueError(f"Alias '{alias}' not found in {fc}")


mean_field = get_alias_not_field_name(test_pts, "MeanValue")
test_arr = arcpy.da.TableToNumPyArray(
    test_pts,
    ["CWS_Points_PWSID", mean_field],
    skip_nulls=True,
)

test_df = pd.DataFrame(test_arr)
test_df.rename(columns={mean_field: "Observed"}, inplace=True)

## initializing the lists for RMSE and MAE

rmse_list = []
mae_list = []

for idx, row in IDW_params.iterrows():

    run_start = time.time()

    raster = row["raster"]
    out_name = row["out_name"]

    zonal_table = os.path.join(fgdb, f"IDW_Zonal_{out_name}")

    ##zonal stats with the 10mi buffers to create mini rasters
    ## creates a zone around each point with mean surrounding values
    ## evens out/ avoids potential sinks in raster when test points extracts values
    try:
        out_ras = arcpy.sa.ZonalStatisticsAsTable(
            in_zone_data=test_bufs5mi,
            zone_field="CWS_Points_PWSID",
            in_value_raster=raster,
            out_table=zonal_table,
            statistics_type="MEAN",
            ignore_nodata="DATA",
        )
    except arcpy.ExecuteError:
        print(f"Run failed: {out_name}")
        print(arcpy.GetMessages(2))
        continue

    zonal_arr = arcpy.da.TableToNumPyArray(
        zonal_table,
        ["CWS_Points_PWSID", "MEAN"],
        skip_nulls=True,
    )

    zonal_df = pd.DataFrame(zonal_arr)
    zonal_df.rename(columns={"MEAN": "Prediction"}, inplace=True)

    merged = test_df.merge(zonal_df, on="CWS_Points_PWSID", how="inner")

    obs = merged["Observed"].to_numpy()
    preds = merged["Prediction"].to_numpy()

    rmse = np.sqrt(np.mean((obs - preds) ** 2))
    mae = np.mean(np.abs(obs - preds))

    rmse_list.append(rmse)
    mae_list.append(mae)

    if arcpy.Exists(zonal_table):
        arcpy.management.Delete(zonal_table)

    run_dur = time.time() - run_start
    print(
        f"[{datetime.now().strftime('%H:%M:%S')}] Run {idx+1} completed in {run_dur:.2f} seconds"
    )


IDW_params["RMSE"] = rmse_list
IDW_params["MAE"] = mae_list


# creating rankings for rmse and MAE and averaging ranks
IDW_runs_rmsesorted = IDW_params.sort_values("RMSE").reset_index(drop=True)
IDW_runs_rmsesorted["RMSE_Rank"] = range(1, len(IDW_runs_rmsesorted) + 1)

IDW_runs_sort = IDW_runs_rmsesorted.sort_values("MAE").reset_index(drop=True)
IDW_runs_sort["MAE_Rank"] = range(1, len(IDW_runs_sort) + 1)

IDW_runs_sort["Ave_Rank"] = (IDW_runs_sort["MAE_Rank"] + IDW_runs_sort["RMSE_Rank"]) / 2
IDW_runs_sorted = IDW_runs_sort.sort_values("Ave_Rank").reset_index(drop=True)


# save to csv
IDW_runs_sorted.to_csv(
    "C:/Duke/Year 2/MP/Interpolation_testing/PFOA_GW/IDW_RankSorted.csv",
    ",",
)


print(
    f"IDW Iterations ran and saved, see final df for the highest average ranking raster"
)


########################################################################################
"""
Universal Kriging Testing - Calculating Stats
"""

########################################################################################

print("starting stat calcs for universal kriging")

mean_field = get_alias_not_field_name(test_pts, "MeanValue")
test_arr = arcpy.da.TableToNumPyArray(
    test_pts,
    ["CWS_Points_PWSID", mean_field],
    skip_nulls=True,
)

test_df = pd.DataFrame(test_arr)
test_df.rename(columns={mean_field: "Observed"}, inplace=True)

## initializing the lists for RMSE and MAE

rmse_list = []
mae_list = []

for idx, row in UK_params.iterrows():

    run_start = time.time()

    raster = row["raster"]
    out_name = row["out_name"]

    zonal_table = os.path.join(fgdb, f"OK_Zonal2_{out_name}")

    ##zonal stats with the 10mi buffers to create mini rasters
    ## creates a zone around each point with mean surrounding values
    ## evens out/ avoids potential sinks in raster when test points extracts values
    try:
        out_ras = arcpy.sa.ZonalStatisticsAsTable(
            in_zone_data=test_bufs5mi,
            zone_field="CWS_Points_PWSID",
            in_value_raster=raster,
            out_table=zonal_table,
            statistics_type="MEAN",
            ignore_nodata="DATA",
        )
    except arcpy.ExecuteError:
        print(f"Run failed: {out_name}")
        print(arcpy.GetMessages(2))
        continue

    zonal_arr = arcpy.da.TableToNumPyArray(
        zonal_table,
        ["CWS_Points_PWSID", "MEAN"],
        skip_nulls=True,
    )

    zonal_df = pd.DataFrame(zonal_arr)
    zonal_df.rename(columns={"MEAN": "Prediction"}, inplace=True)

    merged = test_df.merge(zonal_df, on="CWS_Points_PWSID", how="inner")

    obs = merged["Observed"].to_numpy()
    preds = merged["Prediction"].to_numpy()

    rmse = np.sqrt(np.mean((obs - preds) ** 2))
    mae = np.mean(np.abs(obs - preds))

    rmse_list.append(rmse)
    mae_list.append(mae)

    if arcpy.Exists(zonal_table):
        arcpy.management.Delete(zonal_table)

    run_dur = time.time() - run_start
    print(
        f"[{datetime.now().strftime('%H:%M:%S')}] Run {idx+1} completed in {run_dur:.2f} seconds"
    )


UK_params["RMSE"] = rmse_list
UK_params["MAE"] = mae_list


# creating rankings for rmse and MAE and averaging ranks
UK_runs_rmsesorted = UK_params.sort_values("RMSE").reset_index(drop=True)
UK_runs_rmsesorted["RMSE_Rank"] = range(1, len(UK_runs_rmsesorted) + 1)

UK_runs_sort = UK_runs_rmsesorted.sort_values("MAE").reset_index(drop=True)
UK_runs_sort["MAE_Rank"] = range(1, len(UK_runs_sort) + 1)

UK_runs_sort["Ave_Rank"] = (UK_runs_sort["MAE_Rank"] + UK_runs_sort["RMSE_Rank"]) / 2
UK_runs_sorted = UK_runs_sort.sort_values("Ave_Rank").reset_index(drop=True)


# save to csv
UK_runs_sorted.to_csv(
    "C:/Duke/Year 2/MP/Interpolation_testing/PFOA_GW/UK_RankSorted.csv",
    ",",
)


print(
    f"UK Iterations ran and saved, see final df for the highest average ranking raster"
)


########################################################################################
"""
Ordinary Kriging Testing - Calculating Stats
"""

########################################################################################

print("starting ordinary kriging stat interpretation")

mean_field = get_alias_not_field_name(test_pts, "MeanValue")
test_arr = arcpy.da.TableToNumPyArray(
    test_pts,
    ["CWS_Points_PWSID", mean_field],
    skip_nulls=True,
)

test_df = pd.DataFrame(test_arr)
test_df.rename(columns={mean_field: "Observed"}, inplace=True)

## initializing the lists for RMSE and MAE

rmse_list = []
mae_list = []

for idx, row in OK_params.iterrows():

    run_start = time.time()

    raster = row["raster"]
    out_name = row["out_name"]

    zonal_table = os.path.join(fgdb, f"OK_Zonal_{out_name}")

    ##zonal stats with the 10mi buffers to create mini rasters
    ## creates a zone around each point with mean surrounding values
    ## evens out/ avoids potential sinks in raster when test points extracts values
    try:
        out_ras = arcpy.sa.ZonalStatisticsAsTable(
            in_zone_data=test_bufs5mi,
            zone_field="CWS_Points_PWSID",
            in_value_raster=raster,
            out_table=zonal_table,
            statistics_type="MEAN",
            ignore_nodata="DATA",
        )
    except arcpy.ExecuteError:
        print(f"Run failed: {out_name}")
        print(arcpy.GetMessages(2))
        continue

    zonal_arr = arcpy.da.TableToNumPyArray(
        zonal_table,
        ["CWS_Points_PWSID", "MEAN"],
        skip_nulls=True,
    )

    zonal_df = pd.DataFrame(zonal_arr)
    zonal_df.rename(columns={"MEAN": "Prediction"}, inplace=True)

    merged = test_df.merge(zonal_df, on="CWS_Points_PWSID", how="inner")

    obs = merged["Observed"].to_numpy()
    preds = merged["Prediction"].to_numpy()

    rmse = np.sqrt(np.mean((obs - preds) ** 2))
    mae = np.mean(np.abs(obs - preds))

    rmse_list.append(rmse)
    mae_list.append(mae)

    if arcpy.Exists(zonal_table):
        arcpy.management.Delete(zonal_table)

    run_dur = time.time() - run_start
    print(
        f"[{datetime.now().strftime('%H:%M:%S')}] Run {idx+1} completed in {run_dur:.2f} seconds"
    )


OK_params["RMSE"] = rmse_list
OK_params["MAE"] = mae_list


# creating rankings for rmse and MAE and averaging ranks
OK_runs_rmsesorted = OK_params.sort_values("RMSE").reset_index(drop=True)
OK_runs_rmsesorted["RMSE_Rank"] = range(1, len(OK_runs_rmsesorted) + 1)

OK_runs_sort = OK_runs_rmsesorted.sort_values("MAE").reset_index(drop=True)
OK_runs_sort["MAE_Rank"] = range(1, len(OK_runs_sort) + 1)

OK_runs_sort["Ave_Rank"] = (OK_runs_sort["MAE_Rank"] + OK_runs_sort["RMSE_Rank"]) / 2
OK_runs_sorted = OK_runs_sort.sort_values("Ave_Rank").reset_index(drop=True)


# save to csv
OK_runs_sorted.to_csv(
    "C:/Duke/Year 2/MP/Interpolation_testing/PFOA_GW/OK_RankSorted.csv",
    ",",
)


print(
    f"OK Iterations ran and saved, see final df for the highest average ranking raster"
)

overall_time = time.time() - overall_script_start

print(f"SUCCESS SCRIPT IS DONE in {overall_time/3600:.2f} hours")
