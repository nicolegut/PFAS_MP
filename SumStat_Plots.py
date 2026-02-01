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
# wd_path = Path.cwd().parents[1]
# .chdir(wd_path)

# reading in the RMSE/ MAE data
# subsetting the data

# data = pd.read_csv("./Interpolation_testing/Int_Stats.csv")
# all_sorted_data = data.apply(lambda col: col.sort_values().values)

data = pd.read_csv("./Interpolation_testing/Top15_OrdKrig_MethodComp.csv")

top3_data = data.iloc[0:3, :]


# functions
def long_method_stat_comp(data):
    df_long = data.melt(var_name="method_stat", value_name="value").dropna()
    df_long["method"] = df_long["method_stat"].str.replace(
        r"_(RMSE_ppt|MAE_ppt)", "", regex=True
    )
    df_long["metric"] = df_long["method_stat"].str.extract(r"(RMSE|MAE)")
    return df_long


def plot_interp_RMSE(df, top_df):

    sns.set(style="whitegrid", font_scale=1.2)

    plt.figure(figsize=(8, 6))

    sns.boxplot(
        data=df[df["metric"] == "RMSE"],
        x="method",
        y="value",
        width=0.5,
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
        data=top_df[top_df["metric"] == "RMSE"],
        x="method",
        y="value",
        color="red",
        alpha=1,
        jitter=True,
        size=5,
    )

    plt.ylabel("RMSE (ppt)")
    plt.xlabel("Test Method")
    plt.title("RMSE by Test Method")
    plt.tight_layout()
    return plt.show()


def plot_interp_MAE(df, top_df):

    sns.set(style="whitegrid", font_scale=1.2)

    plt.figure(figsize=(8, 6))

    sns.boxplot(
        data=df[df["metric"] == "MAE"],
        x="method",
        y="value",
        width=0.5,
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
        data=top_df[top_df["metric"] == "MAE"],
        x="method",
        y="value",
        color="red",
        alpha=1,
        jitter=True,
        size=5,
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

data_no_outliers = data.apply(remove_outliers_iqr)


long_out_rm = long_method_stat_comp(data_no_outliers)
plot_interp_MAE(long_out_rm)
plot_interp_RMSE(long_out_rm)


# top 10
top15 = data_no_outliers.iloc[0:15]
long_15 = long_method_stat_comp(top15)
plot_interp_MAE(long_15)
plot_interp_RMSE(long_15)

print(data_no_outliers.notna().sum())


long_data_method_comp = long_method_stat_comp(data)
long_topdata_mcomp = long_method_stat_comp(top3_data)
plot_interp_MAE(long_data_method_comp, long_topdata_mcomp)
plot_interp_RMSE(long_data_method_comp, long_topdata_mcomp)
