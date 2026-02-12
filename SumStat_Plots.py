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

data = pd.read_csv("./Interpolation_testing/PFOA_SW/PFOA_SW_Top15_Comps.csv")

# data = pd.read_csv("./Interpolation_testing/PFOS_SW/PFOS_SW_Top15_Comps.csv")
# data = pd.read_csv("./Interpolation_testing/PFOS_GW/PFOS_GW_Top15_Comps.csv")

data_ukrm = data.iloc[:, :-2]

top5_data = data_ukrm.iloc[0:5, :]
top1_data = data_ukrm.iloc[0:1, :]


# functions
##


def long_method_stat_comp(data):
    df_long = data.melt(var_name="method_stat", value_name="value").dropna()
    df_long["method"] = df_long["method_stat"].str.replace(
        r"_(RMSE_ppt|MAE_ppt)", "", regex=True
    )
    df_long["metric"] = df_long["method_stat"].str.extract(r"(RMSE|MAE)")
    return df_long


def plot_interp_RMSE(df, df2, df3):

    sns.set(style="whitegrid", font_scale=1.3)

    plt.figure(figsize=(8, 6))

    sns.boxplot(
        data=df[df["metric"] == "RMSE"],
        x="method",
        y="value",
        width=0.6,
        showfliers=False,
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
    plt.title("RMSE by Test Method")
    plt.tight_layout()
    return plt.show()


def plot_interp_MAE(df, df2, df3):

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
    plt.title("MAE by Test Method")
    plt.tight_layout()
    return plt.show()


def remove_outliers_iqr(col):
    Q1 = col.quantile(0.25)
    Q3 = col.quantile(0.75)
    IQR = Q3 - Q1
    lower = Q1 - (1.5 * IQR)
    upper = Q3 + (1.5 * IQR)
    return col.where((col >= lower) & (col <= upper))


##outlier removal
# data_no_outliers = data.apply(remove_outliers_iqr)


# for the plots with two levels of coloring
long_15 = long_method_stat_comp(data_ukrm)
long_5 = long_method_stat_comp(top5_data)
long_1 = long_method_stat_comp(top1_data)
plot_interp_MAE(long_15, long_5, long_1)
plot_interp_RMSE(long_15, long_5, long_1)
