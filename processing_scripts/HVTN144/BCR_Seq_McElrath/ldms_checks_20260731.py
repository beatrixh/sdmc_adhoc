import pandas as pd
import numpy as np
import sdmc_tools.access_ldms as access_ldms
import sdmc_tools.process as sdmc
import datetime
import os

df = pd.read_csv(
    "/trials/vaccine/p144/s001/qdata/LabData/BCR_sequencing_pass-through/uploaded_by_lab/McElrath/20260730-04/2026-07-30_HVTN144_BCS2123-BCS2128_AIRR_filtered.tsv",
    sep="\t"
)

# pivot summary
summary= df.groupby(['PTID','Visit']).count().max(axis=1).reset_index()
summary = summary.pivot(
    index='PTID',
    columns='Visit',
    values=0
).fillna(0)

# ldms checks
ldms = access_ldms.pull_one_protocol('hvtn', 144)
ldms['guspec_core'] = ldms.guspec.str.rpartition("-")[0]
ldms = ldms.drop(columns='guspec').drop_duplicates()

df_check = df[['Global_Spec_ID','Visit','PTID']].drop_duplicates().rename(columns={'PTID':'ptid_lab','Global_Spec_ID':'guspec_core','Visit':'visitno_lab'})
df_check.guspec_core = df_check.guspec_core.str.replace("\t","")

df_check = df_check.merge(ldms, on='guspec_core', how='left')
df_check = df_check.rename(columns={'txtpid':'ptid_ldms', 'vidval':'visitno_ldms'})


assert len(df_check.loc[(df_check.visitno_lab.astype(float) != df_check.visitno_ldms.astype(float)),['guspec_core','visitno_lab','visitno_ldms']]) == 0
assert len(df_check.loc[(df_check.ptid_lab.astype(float)!=df_check.ptid_ldms.astype(float)),['guspec_core','ptid_lab','ptid_ldms']]) == 0

# save to .xlsx
savedir = "/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN144/assays/BCR_Sequencing/misc_files/data_processing/"
today = datetime.date.today().isoformat()
summary.to_excel(
    savedir + f'2026-07-30_HVTN144_BCS2123-BCS2128_AIRR_filtered_sample_summary_{today}.xlsx',
    index=True
)
