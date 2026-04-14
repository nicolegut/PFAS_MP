### PFOA GW confusion matrices
## Importing Packages
import pandas as pd
import arcpy, os, time
from datetime import datetime
import numpy as np
from sklearn.metrics import confusion_matrix, accuracy_score, f1_score
import matplotlib.pyplot as plt

# setting up paths
basepath = "C:/Duke/Year 2/MP/"
fgdb = os.path.join(basepath, "PFAS_MP.gdb")
testing_gdb = os.path.join(fgdb, "Testing_Pts")

PFOA_GW_folder = os.path.join(basepath, r"Interpolation_testing\PFOA_GW")
PFOS_GW_folder = os.path.join(basepath, r"Interpolation_testing\PFOS_GW")


# parameter tables
PFOA_GW_IDW = os.path.join(PFOA_GW_folder, "IDW_5_Params.csv")
PFOA_GW_UK = os.path.join(PFOA_GW_folder, "UK_5_Params.csv")
PFOA_GW_OK = os.path.join(PFOA_GW_folder, "OK_5_Params.csv")

PFOS_GW_IDW = os.path.join(PFOS_GW_folder, "IDW_5_Params_sgw.csv")
PFOS_GW_UK = os.path.join(PFOS_GW_folder, "UK_5_Params_sgw.csv")
PFOS_GW_OK = os.path.join(PFOS_GW_folder, "OK_5_Params_sgw.csv")

print("set environments / loading csv")

mapping = {"8000": "8", "16000": "16", "20000": "20", "24000": "24", "#": "Arc"}

# read csvs
IDW_params = pd.read_csv(PFOA_GW_IDW)
IDW_params["out_name"] = IDW_params["raster"].apply(os.path.basename)

OK_params = pd.read_csv(PFOA_GW_OK)
OK_params["out_name"] = OK_params["raster"].apply(os.path.basename)
OK_params["lag_name"] = OK_params["lag_dist"].map(mapping)

UK_params = pd.read_csv(PFOA_GW_UK)
UK_params["out_name"] = UK_params["raster"].apply(os.path.basename)
UK_params["lag_name"] = UK_params["lag_dist"].map(mapping)


# read csvs
IDW_params2 = pd.read_csv(PFOS_GW_IDW)
IDW_params2["out_name"] = IDW_params2["raster"].apply(os.path.basename)

OK_params2 = pd.read_csv(PFOS_GW_OK)
OK_params2["out_name"] = OK_params2["raster"].apply(os.path.basename)
OK_params2["lag_name"] = OK_params2["lag_dist"].map(mapping)

UK_params2 = pd.read_csv(PFOS_GW_UK)
UK_params2["out_name"] = UK_params2["raster"].apply(os.path.basename)
UK_params2["lag_name"] = UK_params2["lag_dist"].map(mapping)


## defining columns
IDW_cols = ["type", "power", "num_points", "distance", "f_pts"]
UK_cols = ["SemiVar_model", "type", "lag_dist", "num_points", "distance", "f_pts"]
OK_cols = ["SemiVar_model", "type", "lag_dist", "num_points", "distance", "f_pts"]

IDW_full = pd.concat([IDW_params, IDW_params2])
OK_full = pd.concat([OK_params, OK_params2])
UK_full = pd.concat([UK_params, UK_params2])

IDW_full = IDW_full[IDW_cols]
OK_full = OK_full[OK_cols]
UK_full = UK_full[UK_cols]


## IDW PARAMS
fig, axes = plt.subplots(2, 3, figsize=(15, 10))

colors = plt.cm.Set2(range(12))

for idx, col in enumerate(IDW_cols):
    ax = axes[idx // 3, idx % 3]
    counts = IDW_full[col].value_counts()
    ax.pie(
        counts.values,
        labels=counts.index,
        autopct="%1.1f%%",
        textprops={"fontsize": 16, "family": "Times New Roman"},
        colors=colors,
    )
    ax.set_title(col, fontsize=20, fontweight="bold", fontfamily="Times New Roman")

axes[1, 2].axis("off")
plt.tight_layout()
plt.show()

##OK PARAMS
fig, axes = plt.subplots(2, 3, figsize=(15, 10))

colors = plt.cm.Set2(range(12))

for idx, col in enumerate(OK_cols):
    ax = axes[idx // 3, idx % 3]
    counts = OK_full[col].value_counts()
    ax.pie(
        counts.values,
        labels=counts.index,
        autopct="%1.1f%%",
        textprops={"fontsize": 16, "family": "Times New Roman"},
        colors=colors,
    )
    ax.set_title(col, fontsize=20, fontweight="bold", fontfamily="Times New Roman")

axes[1, 2].axis("off")
plt.tight_layout()
plt.show()


## UK PARAMS
fig, axes = plt.subplots(2, 3, figsize=(15, 10))

colors = plt.cm.Set2(range(12))

for idx, col in enumerate(UK_cols):
    ax = axes[idx // 3, idx % 3]
    counts = UK_full[col].value_counts()
    ax.pie(
        counts.values,
        labels=counts.index,
        autopct="%1.1f%%",
        textprops={"fontsize": 16, "family": "Times New Roman"},
        colors=colors,
    )
    ax.set_title(col, fontsize=20, fontweight="bold", fontfamily="Times New Roman")

axes[1, 2].axis("off")
plt.tight_layout()
plt.show()
