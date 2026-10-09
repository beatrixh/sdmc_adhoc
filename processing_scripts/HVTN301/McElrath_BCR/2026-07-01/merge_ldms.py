import pandas as pd
import numpy as np
import sdmc_tools.access_ldms as access_ldms
import datetime
import os

input_data_path='/trials/vaccine/p301/s001/qdata/LabData/BCR_sequencing_pass-through/uploaded_by_lab/Hyrien/20260713-03/20260701_HVTN301_airr.tsv'

df = pd.read_csv(
	    input_data_path,
	    sep="\t"
)
df[['PTID','Visit','Global_Spec_ID']].isna().sum()

guspec_path = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN301/assays/BCR_Sequencing/misc_files/McElrath_lab/HVTN301_BCS2082-4_2092_GlobalID_PTID_Visit_CM.xlsx'
guspec_df = pd.read_excel(guspec_path)

guspec_path2 = "/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN301/assays/BCR_Sequencing/misc_files/McElrath_lab/HVTN301_BCRSeq_samples_still_missing_2026-07-14-1057_CM.csv"
guspec_df2 = pd.read_csv(guspec_path2)

guspec_df['guspec_core'] = guspec_df['Global Spec ID'].str.rpartition("-")[0]
guspec_df2['guspec_core'] = guspec_df2['guspec'].str.rpartition("-")[0]

usecols = ['PTID','Visit','guspec_core','Batch']
guspec_df = pd.concat([guspec_df, guspec_df2])[usecols]
guspec_df.groupby(['PTID','Visit']).nunique().max()
guspec_df = guspec_df[['PTID','Visit','guspec_core']].drop_duplicates()

merge_df = df.merge(guspec_df, on=['PTID','Visit'], how='left')


merge_df.loc[merge_df.Global_Spec_ID.isna(), 'Global_Spec_ID'] = merge_df.loc[merge_df.Global_Spec_ID.isna(), 'guspec_core']

today = datetime.date.today().isoformat()
merge_df.drop(columns='guspec_core').to_csv(
    f'/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN301/assays/BCR_Sequencing/misc_files/data_processing/20260701_HVTN301_airr_with_Global_Spec_ID_{today}.txt',
    sep="\t",
    index=False
)
merge_df.drop(columns='guspec_core').to_csv(
    f'/trials/vaccine/p301/s001/qdata/LabData/BCR_sequencing_pass-through/uploaded_by_lab/Hyrien/20260713-03/processed_by_sdmc/20260701_HVTN301_airr_with_Global_Spec_ID_{today}.txt',
    sep="\t",
    index=False
)

input_data_path='/trials/vaccine/p301/s001/qdata/LabData/BCR_sequencing_pass-through/uploaded_by_lab/Hyrien/20260713-03/20260701_HVTN301_wide.tsv'
df = pd.read_csv(
    input_data_path,
    sep="\t",
)

df[['PTID','Visit','Global_Spec_ID']].isna().sum()

merge_df = df.merge(guspec_df, on=['PTID','Visit'], how='left')
merge_df.loc[merge_df.Global_Spec_ID.isna(), 'Global_Spec_ID'] = merge_df.loc[merge_df.Global_Spec_ID.isna(), 'guspec_core']

today = datetime.date.today().isoformat()
merge_df.drop(columns='guspec_core').to_csv(
    f'/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN301/assays/BCR_Sequencing/misc_files/data_processing/20260701_HVTN301_wide_with_Global_Spec_ID_{today}.txt',
    sep="\t",
    index=False
)
merge_df.drop(columns='guspec_core').to_csv(
    f'/trials/vaccine/p301/s001/qdata/LabData/BCR_sequencing_pass-through/uploaded_by_lab/Hyrien/20260713-03/processed_by_sdmc/20260701_HVTN301_wide_with_Global_Spec_ID_{today}.txt',
    sep="\t",
    index=False
)








## create df with missing ids to share to lab -------------------------------------------------------------------------------------------------------------- ##
# df['merge_col'] = df[['PTID','Visit']].astype(str).agg('|'.join, axis=1)
# guspec_df['merge_col'] = guspec_df[['PTID','Visit']].astype(str).agg('|'.join, axis=1)

# df_merge_guspec = df.merge(guspec_df[['guspec','PTID','Visit']], on=['PTID','Visit'], how='left')

# assert df_merge_guspec.shape[0] == df.shape[0]

# df_merge_guspec.loc[(df_merge_guspec.Global_Spec_ID.isna()) & (df_merge_guspec.guspec.notna()),'Global_Spec_ID'] = df_merge_guspec.loc[(df_merge_guspec.Global_Spec_ID.isna()) & (df_merge_guspec.guspec.notna()),'guspec']

# still_missing = df_merge_guspec.loc[df_merge_guspec.Global_Spec_ID.isna(),['PTID','Visit','Batch']].drop_duplicates()

# still_missing.to_csv(
#     '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN301/assays/BCR_Sequencing/misc_files/data_processing/HVTN301_BCRSeq_samples_still_missing_2026-07-14-1057.csv',
#     index=False
# )

# additional_ids = pd.read_excel(
#     '/trials/vaccine/p301/s001/qdata/LabData/BCR_sequencing_pass-through/uploaded_by_lab/McElrath/20260701-04/HVTN301_BCS2082_GlobalID_PTID_Visit.xlsx'
# )

# df_merge_guspec2 = df_merge_guspec.merge(additional_ids, on=['PTID','Visit'], how='left')

# df_merge_guspec2.loc[(df_merge_guspec2.Global_Spec_ID.isna()) & (df_merge_guspec2['Global Spec ID'].notna()),'Global_Spec_ID'] = df_merge_guspec2.loc[(df_merge_guspec2.Global_Spec_ID.isna()) & (df_merge_guspec2['Global Spec ID'].notna()),'Global Spec ID']

# df_merged = df.merge(guspec_df[['guspec','PTID','Visit']], on=['PTID','Visit','Batch','merge_col'], how='left')

# df_merged.loc[(df_merged.Global_Spec_ID.isna()) & (df_merged.guspec.notna()),'Global_Spec_ID'] = df_merged.loc[(df_merged.Global_Spec_ID.isna()) & (df_merged.guspec.notna()),'guspec']

# missings = df_merged.loc[df_merged.Global_Spec_ID.isna(), ['PTID','Visit','Batch']].drop_duplicates()

# missings.to_csv(
#     '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN301/assays/BCR_Sequencing/misc_files/data_processing/HVTN301_BCRSeq_samples_still_missing_2026-07-14.txt',
#     sep="\t",
#     index=False
# )