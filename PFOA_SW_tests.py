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
PFOS_SW_gdb = os.path.join(basepath, "PFOS_SW_Tests.gdb")
PFOS_SW_folder = os.path.join(basepath, r"Interpolation_testing\PFOS_SW")


# paths to points for subsetting
train_pts = os.path.join(testing_gdb, "PFOS_SW_training")
test_pts = os.path.join(testing_gdb, "PFOS_SW_testing")
test_bufs = os.path.join(testing_gdb, "PFOS_SW_test_10mibuf")
test_bufs5mi = os.path.join(testing_gdb, "PFOS_SW_test_5mibuf")
us_mask = os.path.join(fgdb, "ContigUS_Mask")
us_snap_ras = os.path.join(fgdb, "ContigUS_Raster")
SW_PFOS_IDW = os.path.join(PFOS_SW_folder, "PFOS_SW_75IDW.csv")
SW_PFOS_OK = os.path.join(PFOS_SW_folder, "PFOS_SW_75OK.csv")
SW_PFOS_UK = os.path.join(PFOS_SW_folder, "PFOS_SW_75UK.csv")

arcpy.env.overwriteOutput = True
arcpy.env.mask = us_mask
arcpy.env.extent = us_mask
arcpy.env.snapRaster = us_snap_ras
arcpy.env.cellSize = us_snap_ras
# cell size is 5 mi in meters (8046.7m)

print(f"set environments/ loading csv for surface water")

mapping = {"8000": "8", "16000": "16", "20000": "20", "24000": "24", "#": "Arc"}


# read csv for top 75 and edit so that the params match the iterators
IDW_params = pd.read_csv(SW_PFOS_IDW)
IDW_params["out_name"] = IDW_params["raster"].apply(os.path.basename)

OK_params = pd.read_csv(SW_PFOS_OK)
OK_params["out_name"] = OK_params["raster"].apply(os.path.basename)
OK_params["lag_name"] = OK_params["lag_dist"].map(mapping)
OK_params.drop(columns="raster", inplace=True)

UK_params = pd.read_csv(SW_PFOS_UK)
UK_params["out_name"] = UK_params["raster"].apply(os.path.basename)
UK_params["lag_name"] = UK_params["lag_dist"].map(mapping)
UK_params.drop(columns="raster", inplace=True)


# field interpolating
z_field = "MeanValue"

overall_script_start = time.time()


print(
    f"loaded parameter tables/ starting IDW Now [{datetime.now().strftime('%H:%M:%S')}]"
)


########################################################################################
"""
IDW Testing - Generating Rasters
"""
########################################################################################
# list of results to store values before converting to a df
IDW_results = []
IDW_timing = []

IDW_script_start = time.time()

for idx, row in IDW_params.iterrows():

    run_start = time.time()
    out_name = row["out_name"]

    if row["type"] == "variable":
        search_radius = f"VARIABLE {row['num_points']}"
        out_name = out_name

    else:
        search_radius = f"FIXED {row['distance']} {row['f_pts']}"
        out_name = out_name

    # run IDW
    out_ras = arcpy.sa.Idw(
        in_point_features=train_pts,
        z_field=z_field,
        power=row["power"],
        search_radius=search_radius,
    )

    raster_path = os.path.join(PFOS_SW_gdb, out_name)
    out_ras.save(raster_path)

    IDW_results.append(
        [
            raster_path,
            row["out_name"],
            row["type"],
            row["power"],
            row.get("num_points"),
            row.get("distance"),
            row.get("f_pts"),
        ]
    )

    run_dur = time.time() - run_start

    IDW_timing.append(run_dur)

    print(
        f"[{datetime.now().strftime('%H:%M:%S')}] Run {idx+1} completed in {run_dur:.2f} seconds"
    )

total_IDW_dur = time.time() - IDW_script_start
print(f"\nAll {len(IDW_params)} runs completed in {total_IDW_dur/60:.2f} minutes")

## creating a table from the results

IDW_runs = pd.DataFrame(
    IDW_results,
    columns=["raster", "out_name", "type", "power", "num_points", "distance", "f_pts"],
)

IDW_runs.to_csv(
    "C:/Duke/Year 2/MP/Interpolation_testing/PFOS_SW/IDW_t75Raster_Paths.csv", ","
)


########################################################################################
"""
IDW Testing - Calculating Stats
"""
########################################################################################


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

for idx, row in IDW_runs.iterrows():

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


IDW_runs["RMSE"] = rmse_list
IDW_runs["MAE"] = mae_list


# creating rankings for rmse and MAE and averaging ranks
IDW_runs_rmsesorted = IDW_runs.sort_values("RMSE").reset_index(drop=True)
IDW_runs_rmsesorted["RMSE_Rank"] = range(1, len(IDW_runs_rmsesorted) + 1)

IDW_runs_sort = IDW_runs_rmsesorted.sort_values("MAE").reset_index(drop=True)
IDW_runs_sort["MAE_Rank"] = range(1, len(IDW_runs_sort) + 1)

IDW_runs_sort["Ave_Rank"] = (IDW_runs_sort["MAE_Rank"] + IDW_runs_sort["RMSE_Rank"]) / 2
IDW_runs_sorted = IDW_runs_sort.sort_values("Ave_Rank").reset_index(drop=True)


# save to csv
IDW_runs_sorted.to_csv(
    "C:/Duke/Year 2/MP/Interpolation_testing/PFOS_SW/IDW_RankSorted.csv",
    ",",
)


print(
    f"IDW Iterations ran and saved, see final df for the highest average ranking raster"
)


########################################################################################
"""
Univ Krig Testing - Generating Rasters
"""
########################################################################################

UK_results = []
UK_timing = []

UK_script_start = time.time()

for idx, row in UK_params.iterrows():

    out_name = row["out_name"]

    run_start = time.time()

    if row["SemiVar_model"] == "QuadDrift" and row["type"] == "variable":
        kriging_model = f"QuadraticDrift {row['lag_dist']} # # #"
        search_radius = f"VARIABLE {row['num_points']}"
        out_name = f"QuV_lag{row['lag_name']}km_np{str(int(row['num_points']))}"

    elif row["SemiVar_model"] == "LinearDrift" and row["type"] == "variable":
        kriging_model = f"LinearDrift {row['lag_dist']} # # #"
        search_radius = f"VARIABLE {row['num_points']}"
        out_name = f"LnV_lag{row['lag_name']}km_np{str(int(row['num_points']))}"

    elif row["SemiVar_model"] == "QuadDrift" and row["type"] == "fixed":
        kriging_model = f"QuadraticDrift {row['lag_dist']} # # #"
        search_radius = f"FIXED {row['distance']} {row['f_pts']}"
        out_name = f"QuF_lag{row['lag_name']}km_d{str(int(row['distance']/1000))}km_np{str(int(row['f_pts']))}"

    elif row["SemiVar_model"] == "LinearDrift" and row["type"] == "fixed":
        kriging_model = f"LinearDrift {row['lag_dist']} # # #"
        search_radius = f"FIXED {row['distance']} {row['f_pts']}"
        out_name = f"LnF_lag{row['lag_name']}km_d{str(int(row['distance']/1000))}km_np{str(int(row['f_pts']))}"

    else:
        print(f"Error: could not find parameters that matched criteria in row {idx}")
        break

    out_var_raster = os.path.join(PFOS_SW_gdb, f"{out_name}_pdrs")

    # run IDW
    out_ras = arcpy.sa.Kriging(
        in_point_features=train_pts,
        z_field=z_field,
        kriging_model=kriging_model,
        search_radius=search_radius,
        out_variance_prediction_raster=out_var_raster,
    )

    raster_path = os.path.join(PFOS_SW_gdb, out_name)
    out_ras.save(raster_path)

    UK_results.append(
        [
            raster_path,
            row["out_name"],
            row["SemiVar_model"],
            row["type"],
            row["lag_dist"],
            row.get("num_points"),
            row.get("distance"),
            row.get("f_pts"),
        ]
    )

    run_dur = time.time() - run_start
    UK_timing.append(run_dur)

    print(
        f"[{datetime.now().strftime('%H:%M:%S')}] Run {idx+1} completed in {run_dur:.2f} seconds"
    )

total_UK_dur = time.time() - UK_script_start
print(
    f"\nAll {len(UK_params)} Univ Krig raster runs completed in {total_UK_dur/60:.2f} minutes"
)

## creating a table from the results

UK_runs = pd.DataFrame(
    UK_results,
    columns=[
        "raster",
        "out_name",
        "SemiVar_model",
        "type",
        "lag_dist",
        "num_points",
        "distance",
        "f_pts",
    ],
)


##SAVE AS A CSV !!! SO YOU DON'T NEED TO RUN EVERYTHING AGAIN + WAIT 20 MIN
UK_runs.to_csv(
    "C:/Duke/Year 2/MP/Interpolation_testing/PFOS_SW/UK_t75Raster_Paths.csv", ","
)


########################################################################################
"""
Universal Kriging Testing - Calculating Stats
"""

########################################################################################


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

for idx, row in UK_runs.iterrows():

    run_start = time.time()

    raster = row["raster"]
    out_name = row["out_name"]

    zonal_table = os.path.join(fgdb, f"UK_Zonal_{out_name}")

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


UK_runs["RMSE"] = rmse_list
UK_runs["MAE"] = mae_list


# creating rankings for rmse and MAE and averaging ranks
UK_runs_rmsesorted = UK_runs.sort_values("RMSE").reset_index(drop=True)
UK_runs_rmsesorted["RMSE_Rank"] = range(1, len(UK_runs_rmsesorted) + 1)

UK_runs_sort = UK_runs_rmsesorted.sort_values("MAE").reset_index(drop=True)
UK_runs_sort["MAE_Rank"] = range(1, len(UK_runs_sort) + 1)

UK_runs_sort["Ave_Rank"] = (UK_runs_sort["MAE_Rank"] + UK_runs_sort["RMSE_Rank"]) / 2
UK_runs_sorted = UK_runs_sort.sort_values("Ave_Rank").reset_index(drop=True)


# save to csv
UK_runs_sorted.to_csv(
    "C:/Duke/Year 2/MP/Interpolation_testing/PFOS_SW/UK_RankSorted.csv",
    ",",
)


print(
    f"UK Iterations ran and saved, see final df for the highest average ranking raster"
)


########################################################################################
"""
Univ Krig Testing - Generating Rasters
"""
########################################################################################


OK_start = time.time()

# list of results to store values before converting to a df
OK_results = []
OK_timing = []

### IDW Testing

for idx, row in OK_params.iterrows():

    run_start = time.time()
    out_name = row["out_name"]

    ##spherical semivariable fixed vs variable search radii
    if row["SemiVar_model"] == "Spherical" and row["type"] == "variable":
        kriging_model = f"Spherical {row['lag_dist']} # # #"
        search_radius = f"VARIABLE {row['num_points']}"
        out_name = out_name
    elif row["SemiVar_model"] == "Spherical" and row["type"] == "fixed":
        kriging_model = f"Spherical {row['lag_dist']} # # #"
        search_radius = f"FIXED {row['distance']} {row['f_pts']}"
        out_name = out_name

    ##Circular semivariograms fixed vs variable search radii
    elif row["SemiVar_model"] == "Circular" and row["type"] == "variable":
        kriging_model = f"Circular {row['lag_dist']} # # #"
        search_radius = f"VARIABLE {row['num_points']}"
        out_name = out_name
    elif row["SemiVar_model"] == "Circular" and row["type"] == "fixed":
        kriging_model = f"Circular {row['lag_dist']} # # #"
        search_radius = f"FIXED {row['distance']} {row['f_pts']}"
        out_name = out_name

    ##Exponential semivariograms fixed vs variable search radii
    ## skipping the 24k lag dist bc arc crashes --> will need to look into why
    ##adding in the new arc auto
    elif row["SemiVar_model"] == "Exponential" and row["lag_dist"] == "#":
        print("Using ArcGIS default lag distance for Exponential")
        kriging_model = "Exponential # # # #"

        if row["type"] == "variable":
            kriging_model = f"Exponential # # # #"
            search_radius = f"VARIABLE {row['num_points']}"
            out_name = out_name

        elif row["type"] == "fixed":
            kriging_model = f"Exponential # # # #"
            search_radius = f"FIXED {row['distance']} {row['f_pts']}"
            out_name = out_name

    elif row["SemiVar_model"] == "Exponential" and row["lag_dist"] == "24000":
        print(f"Skipping row in Exponential bc of long lag dist")
        continue

    elif row["SemiVar_model"] == "Exponential" and row["type"] == "variable":
        kriging_model = f"Exponential {row['lag_dist']} # # #"
        search_radius = f"VARIABLE {row['num_points']}"
        out_name = out_name
    elif row["SemiVar_model"] == "Exponential" and row["type"] == "fixed":
        kriging_model = f"Exponential {row['lag_dist']} # # #"
        search_radius = f"FIXED {row['distance']} {row['f_pts']}"
        out_name = out_name

    ##Gaussian semivariograms fixed vs variable search radii
    ## skipping the 24k lag dist bc arc crashes --> will need to look into why
    elif row["SemiVar_model"] == "Gaussian" and row["lag_dist"] == "#":
        print("Using ArcGIS default lag distance for Gaussian")

        if row["type"] == "variable":
            kriging_model = f"Gaussian # # # #"
            search_radius = f"VARIABLE {row['num_points']}"
            out_name = out_name

        elif row["type"] == "fixed":
            kriging_model = f"Gaussian # # # #"
            search_radius = f"FIXED {row['distance']} {row['f_pts']}"
            out_name = out_name

    elif row["SemiVar_model"] == "Gaussian" and row["lag_dist"] == "24000":
        print(f"Skipping row in Gaussian bc of long lag dist")
        continue

    elif row["SemiVar_model"] == "Gaussian" and row["type"] == "variable":
        kriging_model = f"Gaussian {row['lag_dist']} # # #"
        search_radius = f"VARIABLE {row['num_points']}"
        out_name = out_name
    elif row["SemiVar_model"] == "Gaussian" and row["type"] == "fixed":
        kriging_model = f"Gaussian {row['lag_dist']} # # #"
        search_radius = f"FIXED {row['distance']} {row['f_pts']}"
        out_name = out_name

    ##Linear semivariograms fixed vs variable search radii
    elif row["SemiVar_model"] == "Linear" and row["type"] == "variable":
        kriging_model = f"Linear {row['lag_dist']} # # #"
        search_radius = f"VARIABLE {row['num_points']}"
        out_name = out_name
    elif row["SemiVar_model"] == "Linear" and row["type"] == "fixed":
        kriging_model = f"Linear {row['lag_dist']} # # #"
        search_radius = f"FIXED {row['distance']} {row['f_pts']}"
        out_name = out_name
    else:
        print(f"Error: could not find parameters that matched criteria in row {idx}")
        break

    out_var_raster = os.path.join(PFOS_SW_gdb, f"{out_name}_pdrs")

    # run kriging
    try:
        out_ras = arcpy.sa.Kriging(
            in_point_features=train_pts,
            z_field=z_field,
            kriging_model=kriging_model,
            search_radius=search_radius,
            out_variance_prediction_raster=out_var_raster,
        )
    except arcpy.ExecuteError:
        print(f"Run failed: {out_name}")
        continue

    raster_path = os.path.join(PFOS_SW_gdb, out_name)
    out_ras.save(raster_path)

    OK_results.append(
        [
            raster_path,
            row["out_name"],
            row["SemiVar_model"],
            row["type"],
            row["lag_dist"],
            row.get("num_points"),
            row.get("distance"),
            row.get("f_pts"),
        ]
    )

    run_dur = time.time() - run_start
    OK_timing.append(run_dur)

    print(
        f"[{datetime.now().strftime('%H:%M:%S')}] Run {idx+1} completed in {run_dur:.2f} seconds"
    )

total_dur = time.time() - OK_start
print(
    f"\nAll {len(OK_params)} Ord Krig Raster runs completed in {total_dur/60:.2f} minutes"
)

## creating a table from the results

OK_runs = pd.DataFrame(
    OK_results,
    columns=[
        "raster",
        "out_name",
        "SemiVar_model",
        "type",
        "lag_dist",
        "num_points",
        "distance",
        "f_pts",
    ],
)


##SAVE AS A CSV !!! SO YOU DON'T NEED TO RUN EVERYTHING AGAIN + WAIT 20 MIN
OK_runs.to_csv(
    "C:/Duke/Year 2/MP/Interpolation_testing/PFOS_SW/OK_t75Raster_Paths.csv", ","
)

########################################################################################
"""
Ordinary Kriging Testing - Calculating Stats
"""

########################################################################################


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

for idx, row in OK_runs.iterrows():

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


OK_runs["RMSE"] = rmse_list
OK_runs["MAE"] = mae_list


# creating rankings for rmse and MAE and averaging ranks
OK_runs_rmsesorted = OK_runs.sort_values("RMSE").reset_index(drop=True)
OK_runs_rmsesorted["RMSE_Rank"] = range(1, len(OK_runs_rmsesorted) + 1)

OK_runs_sort = OK_runs_rmsesorted.sort_values("MAE").reset_index(drop=True)
OK_runs_sort["MAE_Rank"] = range(1, len(OK_runs_sort) + 1)

OK_runs_sort["Ave_Rank"] = (OK_runs_sort["MAE_Rank"] + OK_runs_sort["RMSE_Rank"]) / 2
OK_runs_sorted = OK_runs_sort.sort_values("Ave_Rank").reset_index(drop=True)


# save to csv
OK_runs_sorted.to_csv(
    "C:/Duke/Year 2/MP/Interpolation_testing/PFOS_SW/OK_RankSorted.csv",
    ",",
)


print(
    f"OK Iterations ran and saved, see final df for the highest average ranking raster"
)

overall_time = time.time() - overall_script_start

print(f"SUCCESS SCRIPT IS DONE in {overall_time/3600:.2f} hours")
