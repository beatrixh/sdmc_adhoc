import pandas as pd
import numpy as np
import sdmc_tools.access_ldms as access_ldms
import datetime
import os

input_data_path = '/trials/vaccine/p144/s001/qdata/LabData/BCR_sequencing_pass-through/uploaded_by_lab/McElrath/20260430-01/2026-04-30_HVTN144_BCS2105-BCS2106_AIRR_filtered.tsv'
df = pd.read_csv(input_data_path, sep="\t")


# check that this doesnt have weird things to add/drop
# df.Global_Spec_ID.unique()

# rerun ldms checks
ldms = access_ldms.pull_one_protocol('hvtn', 144)
ldms['guspec_core'] = ldms.guspec.str.rpartition("-")[0]
ldms = ldms.drop(columns='guspec').drop_duplicates()

df_check = df[['Global_Spec_ID','Visit','PTID']].drop_duplicates().rename(columns={'PTID':'ptid_lab','Global_Spec_ID':'guspec_core','Visit':'visitno_lab'})
df_check.guspec_core = df_check.guspec_core.str.replace("\t","")

df_check = df_check.merge(ldms, on='guspec_core', how='left')
df_check = df_check.rename(columns={'txtpid':'ptid_ldms', 'vidval':'visitno_ldms'})

assert len(df_check.loc[(df_check.visitno_lab.astype(float) != df_check.visitno_ldms.astype(float)),['guspec_core','visitno_lab','visitno_ldms']]) == 0
assert len(df_check.loc[(df_check.ptid_lab.astype(float)!=df_check.ptid_ldms.astype(float)),['guspec_core','ptid_lab','ptid_ldms']]) == 0

assert len(df) == df.cell_id.nunique()

summary = pd.pivot_table(
    df,
    index=['PTID'],
    columns='Visit',
    values='cell_id',
    aggfunc='count'
)
savedir = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN144/assays/BCR_Sequencing/misc_files/data_processing/'
today = datetime.date.today().isoformat()
summary.to_excel(
    savedir + f"HVTN144_BCR_Sequencing_McElrath_sample_summary_{today}.xlsx"
)