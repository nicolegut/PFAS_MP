### UCMR4 Data Load
### to see what a final dataset would look like - 4 samples?

import pandas as pd

path4 = "C:/Duke/Year 2/MP/Data/Initial Data/ucmr4_occurrence_data"
UCMR4_All = "UCMR4_All"
UCMR4_data = pd.read_csv(f"{path4}/{UCMR4_All}.txt", sep="\t", encoding="latin1")
