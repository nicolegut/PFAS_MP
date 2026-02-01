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
ordkrig_gdb = os.path.join(basepath, "OrdKrig_Tests2.gdb")
ordkrig_folder = os.path.join(basepath, r"Interpolation_testing\OrdKrig")
# ord_zonal_gdb = os.path.join(basepath, "OrdKrig_ZonalSt.gdb")


# paths to points for subsetting
train_pts = os.path.join(testing_gdb, "PFOA_GW_training")
test_pts = os.path.join(testing_gdb, "PFOA_GW_testing")
test_bufs = os.path.join(testing_gdb, "PFOA_GW_test_10mibuf")
test_bufs5mi = os.path.join(testing_gdb, "PFOA_GW_test_5mibuf")
us_mask = os.path.join(fgdb, "ContigUS_Mask")
us_snap_ras = os.path.join(fgdb, "ContigUS_Raster")
OK_top30_pathscsv = os.path.join(ordkrig_folder, "Ord_top_30.csv")

arcpy.env.overwriteOutput = True
# arcpy.env.mask = us_mask
# arcpy.env.extent = us_mask
# arcpy.env.snapRaster = us_snap_ras
arcpy.env.cellSize = us_snap_ras

print(f"set environments/ loading csv for ord krig")

# read csv instead of reloading all paths
df_runs = pd.read_csv(OK_top30_pathscsv)

df_runs["out_name"] = df_runs["raster"].apply(os.path.basename)


### plan is to create another loop that calculates RMSE from the paths saved in the df_runs df
## would use the path + the testing points path to run extract values from raster
##compare the extracted value (predictions) with the calculated meanvalue (observations)


##RMSE r##uns!!
# calculating rmse for each run/ raster that was created
script_start = time.time()


# ran into errors with alias vs field naming systems -
# I think it will be an issue for both the mean values and predictions so making a function
# https://pro.arcgis.com/en/pro-app/latest/arcpy/functions/gdb-table-properties.htm
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


rmse_list = []
mae_list = []

for idx, row in df_runs.iterrows():

    run_start = time.time()

    raster = row["raster"]
    out_name = row["out_name"]

    ord_zonal_table = os.path.join(fgdb, f"OK_Zonal_{out_name}")

    ##zonal stats with the 10mi buffers to create mini rasters
    ## creates a zone around each point with mean surrounding values
    ## evens out/ avoids potential sinks in raster when test points extracts values
    try:
        out_ras = arcpy.sa.ZonalStatisticsAsTable(
            in_zone_data=test_bufs5mi,
            zone_field="CWS_Points_PWSID",
            in_value_raster=raster,
            out_table=ord_zonal_table,
            statistics_type="MEAN",
            ignore_nodata="DATA",
        )
    except arcpy.ExecuteError:
        print(f"Run failed: {out_name}")
        continue

    zonal_arr = arcpy.da.TableToNumPyArray(
        ord_zonal_table,
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

    run_dur = time.time() - run_start
    print(
        f"[{datetime.now().strftime('%H:%M:%S')}] Run {idx+1} completed in {run_dur:.2f} seconds"
    )


df_runs["RMSE"] = rmse_list
df_runs["MAE"] = mae_list


# creating rankings for rmse and MAE and averaging ranks
df_runs_rmsesorted = df_runs.sort_values("RMSE").reset_index(drop=True)
df_runs_rmsesorted["RMSE_Rank"] = range(1, len(df_runs_rmsesorted) + 1)

df_runs_sort = df_runs_rmsesorted.sort_values("MAE").reset_index(drop=True)
df_runs_sort["MAE_Rank"] = range(1, len(df_runs_sort) + 1)

df_runs_sort["Ave_Rank"] = (df_runs_sort["MAE_Rank"] + df_runs_sort["RMSE_Rank"]) / 2
df_runs_sorted = df_runs_sort.sort_values("Ave_Rank").reset_index(drop=True)


# save to csv
df_runs_sorted.to_csv(
    "C:/Duke/Year 2/MP/Interpolation_testing/OrdKrig/OrdKrig_RasterTests_Rank_v3.csv",
    ",",
)


print(f"Iterations ran and saved, see final df for the highest average ranking raster")
