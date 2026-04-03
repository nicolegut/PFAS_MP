### PFOA SW confusion matrices
## Importing Packages
import pandas as pd
import arcpy, os, time
from datetime import datetime
import numpy as np
from sklearn.metrics import confusion_matrix, accuracy_score, f1_score

# setting up paths
basepath = "C:/Duke/Year 2/MP/"
fgdb = os.path.join(basepath, "PFAS_MP.gdb")
testing_gdb = os.path.join(fgdb, "Testing_Pts")

PFOA_SW_folder = os.path.join(basepath, r"Interpolation_testing\PFOA_SW")

# paths to points for subsetting
test_pts = os.path.join(testing_gdb, "PFOA_SW_testing")
test_bufs5mi = os.path.join(testing_gdb, "PFOA_SW_test_5mibuf")

us_mask = os.path.join(fgdb, "ContigUS_Mask")
us_snap_ras = os.path.join(fgdb, "ContigUS_Raster")

# parameter tables
PFOA_SW_IDW = os.path.join(PFOA_SW_folder, "IDW_5_Params_asw.csv")
PFOA_SW_UK = os.path.join(PFOA_SW_folder, "UK_5_Params_asw.csv")
PFOA_SW_OK = os.path.join(PFOA_SW_folder, "OK_5_Params_asw.csv")

arcpy.env.overwriteOutput = True
arcpy.env.mask = us_mask
arcpy.env.extent = us_mask
arcpy.env.snapRaster = us_snap_ras
arcpy.env.cellSize = us_snap_ras

print("set environments / loading csv")

mapping = {"8000": "8", "16000": "16", "20000": "20", "24000": "24", "#": "Arc"}

# read csvs
IDW_params = pd.read_csv(PFOA_SW_IDW)
IDW_params["out_name"] = IDW_params["raster"].apply(os.path.basename)

OK_params = pd.read_csv(PFOA_SW_OK)
OK_params["out_name"] = OK_params["raster"].apply(os.path.basename)
OK_params["lag_name"] = OK_params["lag_dist"].map(mapping)

UK_params = pd.read_csv(PFOA_SW_UK)
UK_params["out_name"] = UK_params["raster"].apply(os.path.basename)
UK_params["lag_name"] = UK_params["lag_dist"].map(mapping)

overall_script_start = time.time()

########################################################################################
# helper functions
########################################################################################


def get_alias_not_field_name(fc, alias):
    for f in arcpy.ListFields(fc):
        if f.aliasName == alias:
            return f.name
    raise ValueError(f"Alias '{alias}' not found in {fc}")


def classify_val(v):
    if v < 0.004:
        return 0
    elif v <= 0.01:
        return 1
    else:
        return 2


# load observed values
mean_field = get_alias_not_field_name(test_pts, "MeanValue")

test_arr = arcpy.da.TableToNumPyArray(
    test_pts,
    ["CWS_Points_PWSID", mean_field],
    skip_nulls=True,
)

test_df = pd.DataFrame(test_arr)
test_df.rename(columns={mean_field: "Observed"}, inplace=True)

########################################################################################
"""
IDW Confusion Matrices
"""
########################################################################################

acc_list = []
recall_list = []

print(f"Starting IDW confusion matrices [{datetime.now().strftime('%H:%M:%S')}]")

for idx, row in IDW_params.iterrows():

    run_start = time.time()

    raster = row["raster"]
    out_name = row["out_name"]

    zonal_table = os.path.join(fgdb, f"IDW_CM_{out_name}")

    try:
        arcpy.sa.ZonalStatisticsAsTable(
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

    y_true = np.array([classify_val(v) for v in obs])
    y_pred = np.array([classify_val(v) for v in preds])

    cm = confusion_matrix(y_true, y_pred, labels=[0, 1, 2])

    accuracy = accuracy_score(y_true, y_pred)
    recall_contam = (cm[1, 1] + cm[2, 2]) / (cm[1].sum() + cm[2].sum())

    print(f"\n IDW Raster: {out_name}")
    print("Observed (rows) vs Predicted (columns)")
    print("Classes: 0=<0.004, 1=0.004–0.01, 2=>0.01")
    print(cm)

    print(f"Accuracy: {accuracy:.3f}")
    print(f"Contaminanted Recall: {recall_contam:.3f}")

    acc_list.append(accuracy)
    recall_list.append(recall_contam)

    if arcpy.Exists(zonal_table):
        arcpy.management.Delete(zonal_table)

    run_dur = time.time() - run_start
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Completed in {run_dur:.2f} sec")


mean_acc = np.nanmean(acc_list)
mean_recall = np.nanmean(recall_list)

print("\n=================================")
print("IDW Average Performance")
print(f"Mean Accuracy: {mean_acc:.3f}")
print(f"Mean Recall (>0.004): {mean_recall:.3f}")
print("=================================\n")

########################################################################################
"""
Universal Kriging Confusion Matrices
"""
########################################################################################

acc_list = []
recall_list = []

print(
    f"\nStarting Universal Kriging confusion matrices [{datetime.now().strftime('%H:%M:%S')}]"
)

for idx, row in UK_params.iterrows():

    run_start = time.time()

    raster = row["raster"]
    out_name = row["out_name"]

    zonal_table = os.path.join(fgdb, f"UK_CM_{out_name}")

    try:
        arcpy.sa.ZonalStatisticsAsTable(
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

    y_true = np.array([classify_val(v) for v in obs])
    y_pred = np.array([classify_val(v) for v in preds])

    cm = confusion_matrix(y_true, y_pred, labels=[0, 1, 2])

    accuracy = accuracy_score(y_true, y_pred)
    recall_contam = (cm[1, 1] + cm[2, 2]) / (cm[1].sum() + cm[2].sum())

    print(f"\n UK Raster: {out_name}")
    print("Observed (rows) vs Predicted (columns)")
    print("Classes: 0=<0.004, 1=0.004–0.01, 2=>0.01")
    print(cm)

    print(f"Accuracy: {accuracy:.3f}")
    print(f"Contaminanted Recall: {recall_contam:.3f}")

    acc_list.append(accuracy)
    recall_list.append(recall_contam)

    if arcpy.Exists(zonal_table):
        arcpy.management.Delete(zonal_table)

    run_dur = time.time() - run_start
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Completed in {run_dur:.2f} sec")


mean_acc = np.nanmean(acc_list)
mean_recall = np.nanmean(recall_list)

print("\n=================================")
print("UK Average Performance")
print(f"Mean Accuracy: {mean_acc:.3f}")
print(f"Mean Recall (>0.004): {mean_recall:.3f}")
print("=================================\n")

########################################################################################
"""
Ordinary Kriging Confusion Matrices
"""
########################################################################################

acc_list = []
recall_list = []

print(
    f"\nStarting Ordinary Kriging confusion matrices [{datetime.now().strftime('%H:%M:%S')}]"
)

for idx, row in OK_params.iterrows():

    run_start = time.time()

    raster = row["raster"]
    out_name = row["out_name"]

    zonal_table = os.path.join(fgdb, f"OK_CM_{out_name}")

    try:
        arcpy.sa.ZonalStatisticsAsTable(
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

    y_true = np.array([classify_val(v) for v in obs])
    y_pred = np.array([classify_val(v) for v in preds])

    cm = confusion_matrix(y_true, y_pred, labels=[0, 1, 2])

    accuracy = accuracy_score(y_true, y_pred)
    recall_contam = (cm[1, 1] + cm[2, 2]) / (cm[1].sum() + cm[2].sum())

    print(f"\n OK Raster: {out_name}")
    print("Observed (rows) vs Predicted (columns)")
    print("Classes: 0=<0.004, 1=0.004–0.01, 2=>0.01")
    print(cm)

    print(f"Accuracy: {accuracy:.3f}")
    print(f"Contaminanted Recall: {recall_contam:.3f}")

    acc_list.append(accuracy)
    recall_list.append(recall_contam)

    if arcpy.Exists(zonal_table):
        arcpy.management.Delete(zonal_table)

    run_dur = time.time() - run_start
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Completed in {run_dur:.2f} sec")

overall_time = time.time() - overall_script_start
print(f"\nSUCCESS SCRIPT FINISHED in {overall_time/3600:.2f} hours")


mean_acc = np.nanmean(acc_list)
mean_recall = np.nanmean(recall_list)

print("\n=================================")
print("OK Average Performance")
print(f"Mean Accuracy: {mean_acc:.3f}")
print(f"Mean Recall (>0.004): {mean_recall:.3f}")
print("=================================\n")
