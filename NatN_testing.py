###this is the script for an initial testing of parameters for Natural Neighbor
## use this as the submission for

##Importing Packages
import pandas as pd
import arcpy, os, time
from datetime import datetime
import numpy as np


# setting up paths
basepath = "C:/Duke/Year 2/MP/"
fgdb = os.path.join(basepath, "PFAS_MP.gdb")
contig_gdb = os.path.join(fgdb, "Exp_Int_pts")
testing_gdb = os.path.join(fgdb, "Testing_Pts")
# change this one for different interpolations
natn_gdb = os.path.join(basepath, "NatN_Tests.gdb")
natn_folder = os.path.join(basepath, "Interpolation_testing\\NatN")

# paths to points for subsetting
train_pts = os.path.join(testing_gdb, "PFOA_GW_training")
test_pts = os.path.join(testing_gdb, "PFOA_GW_testing")
us_mask = os.path.join(fgdb, "ContigUS_Mask")

arcpy.env.overwriteOutput = True
arcpy.env.mask = us_mask
arcpy.env.extent = us_mask
arcpy.env.outputCoordinateSystem = arcpy.Describe(train_pts).spatialReference

# functions


def get_alias_not_field_name(fc, alias):
    for f in arcpy.ListFields(fc):
        if f.aliasName == alias:
            return f.name
    raise ValueError(f"Alias '{alias}' not found in {fc}")


# natural neighbor interpolation

z_field = "MeanValue"
cell_size = 8046.7
out_name = "NatN_Int"


script_start = time.time()

# list of results to store values before converting to a df
results = []

out_ras = arcpy.sa.NaturalNeighbor(
    in_point_features=train_pts, z_field=z_field, cell_size=cell_size
)

raster_path = os.path.join(natn_gdb, out_name)
out_ras.save(raster_path)

results.append([raster_path])


total_dur = time.time() - script_start
print(
    f"[{datetime.now().strftime('%H:%M:%S')}] Run completed in {total_dur:.2f} seconds"
)

## creating a table from the results

df_runs = pd.DataFrame(
    results,
    columns=["raster"],
)

rmse_list = []
mae_list = []

field_names = [f.name for f in arcpy.ListFields(test_pts)]

if "Prediction" in field_names:
    arcpy.management.DeleteField(test_pts, "Prediction")

# extract predictions to test points
arcpy.sa.ExtractMultiValuesToPoints(
    in_point_features=test_pts,
    in_rasters=[[raster_path, "Prediction"]],
    bilinear_interpolate_values="NONE",
)

preds = []
obs = []

mean_field = get_alias_not_field_name(test_pts, "MeanValue")
pred_field = get_alias_not_field_name(test_pts, "Prediction")

with arcpy.da.SearchCursor(test_pts, [mean_field, pred_field]) as cur:
    for o, p in cur:
        if p is not None:
            preds.append(p)
            obs.append(o)

rmse = np.sqrt(sum(((np.array(obs) - np.array(preds)) ** 2)) / len(obs))
mae = np.mean(np.abs(np.array(obs) - np.array(preds)))

print(f"{rmse} {raster_path}")
rmse_list.append(rmse)

print(f"{mae} {raster_path}")
mae_list.append(mae)

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
    "C:/Duke/Year 2/MP/Interpolation_testing/NatN/NatN_AveRank_Sort.csv", ","
)


print(
    f"The Natural Neighbor interpolation has an RMSE of {df_runs_sorted['RMSE'].iloc[0]} and a MAE of {df_runs_sorted['MAE'].iloc[0]}"
)
