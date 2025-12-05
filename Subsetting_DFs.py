##subsetting the points for training/ testing


##Importing Packages
import pandas as pd
import arcpy, os, time


# setting up paths
basepath = "C:/Duke/Year 2/MP/"
fgdb = os.path.join(basepath, "PFAS_MP.gdb")
contig_gdb = os.path.join(fgdb, "Exp_Int_pts")
testing_gdb = os.path.join(fgdb, "Testing_Pts")

# paths to points for subsetting
PFOA_GW = os.path.join(contig_gdb, "PFOA_GW_Dec2025_Contig")
PFOA_SW = os.path.join(contig_gdb, "PFOA_SW_Dec2025_Contig")
PFOS_GW = os.path.join(contig_gdb, "PFOS_GW_Dec2025_Contig")
PFOS_SW = os.path.join(contig_gdb, "PFOS_SW_Dec2025_Contig")

point_list = [PFOA_GW, PFOS_GW, PFOA_SW, PFOS_SW]

naming_list = ["PFOA_GW", "PFOS_GW", "PFOA_SW", "PFOS_SW"]

arcpy.env.overwriteOutput = True

for point, name in zip(point_list, naming_list):
    out_train = os.path.join(testing_gdb, f"{name}_training")
    out_test = os.path.join(testing_gdb, f"{name}_testing")

    arcpy.SubsetFeatures_ga(
        in_features=f"{point}",
        out_training_feature_class=out_train,
        out_test_feature_class=out_test,
        size_of_training_dataset="90",
        subset_size_units="PERCENTAGE_OF_INPUT",
    )
