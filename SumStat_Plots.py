## Distributions of Interpolation RMSE and MAE values

# importing libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import os
import seaborn as sns

# setting working directory
# moving back 2 folders to access the MP folder
# current workspace: 'c:/Duke/Year 2/MP/Scripts/PFAS_MP'
# new workspace: 'c:/Duke/Year 2/MP'
"""
Run these two lines below to set the current working directory!! 
Then comment back out to avoid pushing the wd back too much
"""

# wd_path = Path.cwd().parents[1]
# os.chdir(wd_path)

# reading in the RMSE/ MAE data
# subsetting the data

# data = pd.read_csv("./Interpolation_testing/Int_Stats.csv")
# all_sorted_data = data.apply(lambda col: col.sort_values().values)

PFOA_SW_data = pd.read_csv("./Interpolation_testing/PFOA_SW/PFOA_SW_Top15_Comps.csv")
PFOS_SW_data = pd.read_csv("./Interpolation_testing/PFOS_SW/PFOS_SW_Top15_Comps.csv")
PFOS_GW_data = pd.read_csv("./Interpolation_testing/PFOS_GW/PFOS_GW_Top15_Comps.csv")
PFOA_GW_data = pd.read_csv("./Interpolation_testing/PFOA_GW/PFOA_GW_Top15_Comps.csv")

PFOA_SW_nouk = PFOA_SW_data.iloc[:, :-2]
PFOS_SW_nouk = PFOS_SW_data.iloc[:, :-2]
PFOS_GW_nouk = PFOS_GW_data.iloc[:, :-2]
PFOA_GW_nouk = PFOA_GW_data.iloc[:, :-2]

color_palette = sns.color_palette("Set2", 4)


def top_5(data):
    top5 = data.iloc[0:5, :]
    return top5


def top_1(data):
    top5 = data.iloc[0:1, :]
    return top5


def long_method_stat_comp(data):
    df_long = data.melt(var_name="method_stat", value_name="value").dropna()
    df_long["method"] = df_long["method_stat"].str.replace(
        r"_(RMSE_ppt|MAE_ppt)", "", regex=True
    )
    df_long["metric"] = df_long["method_stat"].str.extract(r"(RMSE|MAE)")
    return df_long


def plot_interp_RMSE(df, df2, df3, title):

    sns.set(style="whitegrid", font_scale=1.3)

    plt.figure(figsize=(8, 6))

    sns.boxplot(
        data=df[df["metric"] == "RMSE"],
        x="method",
        y="value",
        width=0.6,
        showfliers=False,
        color="green",
    )

    sns.stripplot(
        data=df[df["metric"] == "RMSE"],
        x="method",
        y="value",
        color="black",
        alpha=0.6,
        jitter=True,
        size=5,
    )
    # top 3
    sns.stripplot(
        data=df2[df2["metric"] == "RMSE"],
        x="method",
        y="value",
        color="orange",
        alpha=1,
        jitter=True,
        size=7,
    )

    sns.stripplot(
        data=df3[df3["metric"] == "RMSE"],
        x="method",
        y="value",
        color="red",
        alpha=1,
        jitter=True,
        size=9,
    )

    plt.ylabel("RMSE (ppt)")
    plt.xlabel("Test Method")
    plt.title(f"RMSE by Test Method: {title}")
    plt.tight_layout()
    return plt.show()


def plot_interp_MAE(df, df2, df3, title):

    sns.set(style="whitegrid", font_scale=1.3)

    plt.figure(figsize=(8, 6))

    sns.boxplot(
        data=df[df["metric"] == "MAE"],
        x="method",
        y="value",
        width=0.6,
        showfliers=False,
    )

    sns.stripplot(
        data=df[df["metric"] == "MAE"],
        x="method",
        y="value",
        color="black",
        alpha=0.6,
        jitter=True,
        size=5,
    )

    sns.stripplot(
        data=df2[df2["metric"] == "MAE"],
        x="method",
        y="value",
        color="orange",
        alpha=1,
        jitter=True,
        size=7,
    )

    sns.stripplot(
        data=df3[df3["metric"] == "MAE"],
        x="method",
        y="value",
        color="red",
        alpha=1,
        jitter=True,
        size=9,
    )

    plt.ylabel("MAE (ppt)")
    plt.xlabel("Test Method")
    plt.title(f"MAE by Test Method: {title}")
    plt.tight_layout()
    return plt.show()


##PFOA GW RMSE/ MAE
PFOA_GW5 = top_5(PFOA_GW_data)
PFOA_GW1 = top_1(PFOA_GW_data)
PFOA_GW_long = long_method_stat_comp(PFOA_GW_data)
PFOA_GW5_long = long_method_stat_comp(PFOA_GW5)
PFOA_GW1_long = long_method_stat_comp(PFOA_GW1)

plot_interp_MAE(PFOA_GW_long, PFOA_GW5_long, PFOA_GW1_long, "PFOA GW")
plot_interp_RMSE(PFOA_GW_long, PFOA_GW5_long, PFOA_GW1_long, "PFOA GW")

##PFOA SW RMSE/ MAE
PFOA_SW5 = top_5(PFOA_SW_data)
PFOA_SW1 = top_1(PFOA_SW_data)
PFOA_SW_long = long_method_stat_comp(PFOA_SW_data)
PFOA_SW5_long = long_method_stat_comp(PFOA_SW5)
PFOA_SW1_long = long_method_stat_comp(PFOA_SW1)

plot_interp_MAE(PFOA_SW_long, PFOA_SW5_long, PFOA_SW1_long, "PFOA SW")
plot_interp_RMSE(PFOA_SW_long, PFOA_SW5_long, PFOA_SW1_long, "PFOA SW")


##PFOS SW RMSE/ MAE
PFOS_SW5 = top_5(PFOS_SW_data)
PFOS_SW1 = top_1(PFOS_SW_data)
PFOS_SW_long = long_method_stat_comp(PFOS_SW_data)
PFOS_SW5_long = long_method_stat_comp(PFOS_SW5)
PFOS_SW1_long = long_method_stat_comp(PFOS_SW1)

plot_interp_MAE(PFOS_SW_long, PFOS_SW5_long, PFOS_SW1_long, "PFOS SW")
plot_interp_RMSE(PFOS_SW_long, PFOS_SW5_long, PFOS_SW1_long, "PFOS SW")

##PFOS GW RMSE/ MAE
PFOS_GW5 = top_5(PFOS_GW_data)
PFOS_GW1 = top_1(PFOS_GW_data)
PFOS_GW_long = long_method_stat_comp(PFOS_GW_data)
PFOS_GW5_long = long_method_stat_comp(PFOS_GW5)
PFOS_GW1_long = long_method_stat_comp(PFOS_GW1)

plot_interp_MAE(PFOS_GW_long, PFOS_GW5_long, PFOS_GW1_long, "PFOS GW")
plot_interp_RMSE(PFOS_GW_long, PFOS_GW5_long, PFOS_GW1_long, "PFOS GW")


PFOS_GW_Timing = pd.read_csv("./Interpolation_testing/PFOS_GW/PFOS_GW_timing_long.csv")
PFOS_GW_Timing = PFOS_GW_Timing.iloc[:, 1:3]

PFOS_SW_Timing = pd.read_csv("./Interpolation_testing/PFOS_SW/Methods_timing_long.csv")


PFOS_GW_Timing["Method"] = (
    PFOS_GW_Timing["variable"]
    .str.replace("_timing", "", regex=False)
    .replace({"IDW": "IDW", "OK": "OrdKrig", "UK": "UnivKrig"})
)

PFOS_SW_Timing["Method"] = PFOS_SW_Timing["Method"].replace(
    {"IDW": "IDW", "Ord_Krig": "OrdKrig", "Univ_Krig": "UnivKrig"}
)

PFOS_GW_Timing = PFOS_GW_Timing.rename(columns={"value": "Time"})
PFOS_SW_Timing = PFOS_SW_Timing.rename(columns={"Time (s)": "Time"})

PFOS_GW_Timing["Source"] = "GW"
PFOS_SW_Timing["Source"] = "SW"

gw = PFOS_GW_Timing[["Method", "Time", "Source"]]
sw = PFOS_SW_Timing[["Method", "Time", "Source"]]

PFOS_Timing_All = pd.concat([gw, sw], ignore_index=True)

plt.figure(figsize=(8, 6))

# ONE box per method (combined GW + SW)
sns.boxplot(
    data=PFOS_Timing_All, x="Method", y="Time", showfliers=False, color="lightgray"
)

# Colored dots by source
sns.stripplot(
    data=PFOS_Timing_All,
    x="Method",
    y="Time",
    hue="Source",
    dodge=False,
    jitter=True,
    alpha=0.7,
    size=6,
)

plt.ylabel("Time (s)")
plt.title("Interpolation Runtime by Method")

plt.legend(title="Water Type")
plt.tight_layout()
plt.show()
