###this is the script for an initial testing of parameters for IDW
## use this as the submission for

##Importing Packages
import pandas as pd
import arcpy, os, time
from datetime import datetime
import itertools
import numpy as np


# setting up paths
basepath = "C:/Duke/Year 2/MP/"
fgdb = os.path.join(basepath, "PFAS_MP.gdb")
contig_gdb = os.path.join(fgdb, "Exp_Int_pts")
testing_gdb = os.path.join(fgdb, "Testing_Pts")
# change this one for different interpolations
idw_gdb = os.path.join(basepath, "IDW_Tests.gdb")

# paths to points for subsetting
train_pts = os.path.join(testing_gdb, "PFOA_GW_training")
test_pts = os.path.join(testing_gdb, "PFOA_GW_testing")
us_mask = os.path.join(fgdb, "ContigUS_Mask")

arcpy.env.overwriteOutput = True
arcpy.env.mask = us_mask
arcpy.env.extent = us_mask

# field interpolating
z_field = "MeanValue"
# power input - general range is 0-3, default is 2
power_inp = [0.5, 1, 2, 3]
#### variable search radius - RadiusVariable ({numberofPoints}, {maxDistance})
var_rad_pts = [5, 12, 50, 100]
# will be using the default for max distance

#### fixed search radius - RadiusFixed ({distance}, {minNumberofPoints})
f_dist = [16000, 40250, 80500, 161000]
# default is 5x cell size of output (5 mi -> 25 mi), going to use 10 / 25/ 50 / 100 mi -
# converted to meter for map units
f_pts = [2, 6, 12, 24]
##tried with the default 0 at first, but needs at least 2 points to fully interpolate the map
## as of 12/5/25 tests with 24 points have not been run

cell_size = 8046.7  # 5 miles in meters - previously used 10 miles, trying this out?

##Ranges of parameters are in lists above
## creating tables of variable vs fixed search distance parameter combinations
## will eventually be a table of all possible combinations/ settings
## loop will iterate through the df for input parameter values


# creating a table of variable radius combinations
var_combos = list(itertools.product(power_inp, var_rad_pts))

df_var = pd.DataFrame(var_combos, columns=["power", "num_points"])
df_var["type"] = "variable"

# creating a table of fixed radius combinations
fixed_combos = list(itertools.product(power_inp, f_dist, f_pts))

df_fix = pd.DataFrame(fixed_combos, columns=["power", "rad_dist", "f_pts"])
df_fix["type"] = "fixed"

# final table of input combinations
idw_combos = pd.concat([df_fix, df_var], axis=0, ignore_index=True)

# time stamps for how long each iteration/ whole loop takes to run
script_start = time.time()

# list of results to store values before converting to a df
results = []


### IDW Testing

for idx, row in idw_combos.iterrows():

    run_start = time.time()

    if row["type"] == "variable":
        search_radius = f"VARIABLE {row['num_points']}"
        out_name = f"V_p{str(int(row['power'])).replace('.','')}_np{str(int(row['num_points']))}"

    else:
        search_radius = f"FIXED {row['rad_dist']} {row['f_pts']}"
        out_name = f"F_p{str(int(row['power'])).replace('.','')}_d{str(int(row['rad_dist']/1000)).replace('.','')}km_np{str(int(row['f_pts']))}"

    # run IDW
    out_ras = arcpy.sa.Idw(
        in_point_features=train_pts,
        z_field=z_field,
        cell_size=cell_size,
        power=row["power"],
        search_radius=search_radius,
    )

    raster_path = os.path.join(idw_gdb, out_name)
    out_ras.save(raster_path)

    results.append(
        [
            raster_path,
            row["type"],
            row["power"],
            row.get("num_points"),
            row.get("rad_dist"),
            row.get("f_pts"),
        ]
    )

    run_dur = time.time() - run_start
    print(
        f"[{datetime.now().strftime('%H:%M:%S')}] Run {idx+1} completed in {run_dur:.2f} seconds"
    )

total_dur = time.time() - script_start
print(f"\nAll {len(idw_combos)} runs completed in {total_dur/60:.2f} minutes")

## creating a table from the results

df_runs = pd.DataFrame(
    results,
    columns=["raster", "type", "power", "num_points", "distance", "f_pts"],
)

### plan is to create another loop that calculates RMSE from the paths saved in the df_runs df
## would use the path + the testing points path to run extract values from raster
##compare the extracted value (predictions) with the calculated meanvalue (observations)

##RMSE r##uns!!
# calculating rmse for each run/ raster that was created
rmse_list = []


for i, row in df_runs.iterrows():
    raster = row["raster"]

    # extract predictions to test points
    arcpy.sa.ExtractMultiValuesToPoints(
        in_point_features=test_pts,
        in_rasters=f"{raster} Prediction",
        bilinear_interpolate_values="NONE",
    )

    preds = []
    obs = []

    with arcpy.da.SearchCursor(test_pts, ["MeanValue", "Pred"]) as cur:
        for o, p in cur:
            if p is not None:
                preds.append(p)
                obs.append(o)

    rmse = np.sqrt(((np.array(obs) - np.array(preds)) ** 2).mean())
    rmse_list.append(rmse)

df_runs["RMSE"] = rmse_list
