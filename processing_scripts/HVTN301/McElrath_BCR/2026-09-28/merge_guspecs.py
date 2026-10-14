import pandas as pd
import numpy as np
import sdmc_tools.access_ldms as access_ldms
import datetime
import os

# read in ju's outputs -------------------------------------------- #
input_data_path='/trials/vaccine/p301/s001/qdata/LabData/BCR_sequencing_pass-through/uploaded_by_lab/Hyrien/20261007-03/20260928_HVTN301_airr.tsv'

df = pd.read_csv(
	    input_data_path,
	    sep="\t"
)
# df[['PTID','Visit','Global_Spec_ID']].isna().sum()



# read in guspec maps -------------------------------------------- #
guspec_path = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN301/assays/BCR_Sequencing/misc_files/McElrath_lab/HVTN301_BCS2082-4_2092_GlobalID_PTID_Visit_CM.xlsx'
guspec_df = pd.read_excel(guspec_path)

guspec_path2 = "/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN301/assays/BCR_Sequencing/misc_files/McElrath_lab/HVTN301_BCRSeq_samples_still_missing_2026-07-14-1057_CM.csv"
guspec_df2 = pd.read_csv(guspec_path2)

guspec_df['guspec_core'] = guspec_df['Global Spec ID'].str.rpartition("-")[0]
guspec_df2['guspec_core'] = guspec_df2['guspec'].str.rpartition("-")[0]

usecols = ['PTID','Visit','guspec_core','Batch']
guspec_df = pd.concat([guspec_df, guspec_df2])[usecols]
guspec_df.Visit = guspec_df.Visit.astype(float)

guspec_df = guspec_df[['PTID','Visit','guspec_core']].drop_duplicates()

# merge on guspec_cores  -------------------------------------------- #
merge_df = df.merge(guspec_df, on=['PTID','Visit'], how='left')
merge_df.loc[merge_df.Global_Spec_ID.isna(), 'Global_Spec_ID'] = merge_df.loc[merge_df.Global_Spec_ID.isna(), 'guspec_core']
assert merge_df.Global_Spec_ID.isna().sum() == 0

merge_df = merge_df.drop(columns='guspec_core')


# save to txt ---------------------------------------------------- #
savedir = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN301/assays/BCR_Sequencing/misc_files/data_processing/'
merge_df.to_csv(
    savedir + '20260928_HVTN301_airr_with_Global_Spec_ID_2026-10-14.txt',
    sep="\t",
    index=False
)

## do the same on wide ---------------------------------------------------- #
input_data_path='/trials/vaccine/p301/s001/qdata/LabData/BCR_sequencing_pass-through/uploaded_by_lab/Hyrien/20261007-03/20260928_HVTN301_wide.tsv'
wide = pd.read_csv(
	    input_data_path,
	    sep="\t"
)

wide_merge = wide.merge(guspec_df, on=['PTID','Visit'], how='left')
wide_merge.loc[wide_merge.Global_Spec_ID.isna(), 'Global_Spec_ID'] = wide_merge.loc[wide_merge.Global_Spec_ID.isna(), 'guspec_core']
assert wide_merge.Global_Spec_ID.isna().sum() == 0
wide_merge = wide_merge.drop(columns='guspec_core')

wide_merge.to_csv(
    savedir + '20260928_HVTN301_wide_with_Global_Spec_ID_2026-10-14.txt',
    sep="\t",
    index=False
)