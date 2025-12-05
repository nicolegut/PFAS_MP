# this script compares the distributions of contaminants/ point values
# to see if only using one df to test them all is okay
# use miniforge env to run this one!


# importing packages
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import scikit_posthocs as sp

# path to the csvs
path = "C:/Duke/Year 2/MP/Data/stats_csvs/"

PFOA_GW = pd.read_csv(f"{path}PFOA_GW_stats.csv")
PFOA_SW = pd.read_csv(f"{path}PFOA_SW_stats.csv")
PFOS_GW = pd.read_csv(f"{path}PFOS_GW_stats.csv")
PFOS_SW = pd.read_csv(f"{path}PFOS_SW_stats.csv")

# kruskall wallace for the means of the groups - using for the testing
# kruskall vs anova bc non-normally distrib
f_stat2, p_val2 = stats.kruskal(
    PFOA_GW["MeanValue"],
    PFOA_SW["MeanValue"],
    PFOS_GW["MeanValue"],
    PFOS_SW["MeanValue"],
)

print(f"the f statistic is {f_stat2} and the p-value is {p_val2}")
# f stat 25.73 and p-val is 1.01*e-5


# tukey to see what the issues/ differences are
# converting to a single long df for the test
all_values = pd.concat(
    [
        PFOA_GW["MeanValue"],
        PFOA_SW["MeanValue"],
        PFOS_GW["MeanValue"],
        PFOS_SW["MeanValue"],
    ],
    ignore_index=True,
)

groups = (
    ["PFOA_GW"] * len(PFOA_GW)
    + ["PFOA_SW"] * len(PFOA_SW)
    + ["PFOS_GW"] * len(PFOS_GW)
    + ["PFOS_SW"] * len(PFOS_SW)
)

# converting to series for the concat
all_values = pd.Series(all_values)
groups = pd.Series(groups)

# concat + renaming the columns
dunn_test = pd.concat([groups, all_values], axis=1)
dunn_test.columns = ["group", "value"]


# tukey test using the long df
dunn = sp.posthoc_dunn(dunn_test, val_col="value", group_col="group")

# creating a boolean for rejecting null hyp
dunn_bool_clean = dunn < 0.05
