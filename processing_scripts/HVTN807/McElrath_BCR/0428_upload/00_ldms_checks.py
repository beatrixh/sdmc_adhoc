import pandas as pd
import numpy as np
import sdmc_tools.access_ldms as access_ldms
import sdmc_tools.process as sdmc
import datetime
import os

input_data_path = '/trials/vaccine/p807/s001/qdata/LabData/BCR_sequencing_pass-through/uploaded_by_lab/McElrath/20260428-01/2026-04-28_HVTN807_BCS2107-BCS2108_AIRR_filtered.tsv'
df = pd.read_csv(input_data_path, sep="\t")

# pivot summary
summary= df.groupby(['PTID','Visit']).count().max(axis=1).reset_index()
summary = summary.pivot(
    index='PTID',
    columns='Visit',
    values=0
).fillna(0)

# ldms checks
ldms = access_ldms.pull_one_protocol('hvtn', 807)
ldms['guspec_core'] = ldms.guspec.str.rpartition("-")[0]
ldms = ldms.drop(columns='guspec').drop_duplicates()

df_check = df[['Global_Spec_ID','Visit','PTID']].drop_duplicates().rename(columns={'PTID':'ptid_lab','Global_Spec_ID':'guspec_core','Visit':'visitno_lab'})
df_check.guspec_core = df_check.guspec_core.str.replace("\t","")

df_check = df_check.merge(ldms, on='guspec_core', how='left')
df_check = df_check.rename(columns={'txtpid':'ptid_ldms', 'vidval':'visitno_ldms'})

assert len(df_check.loc[(df_check.visitno_lab.astype(float) != df_check.visitno_ldms.astype(float)),['guspec_core','visitno_lab','visitno_ldms']]) == 0
assert len(df_check.loc[(df_check.ptid_lab.astype(float)!=df_check.ptid_ldms.astype(float)),['guspec_core','ptid_lab','ptid_ldms']]) == 0

# save to .xlsx
today = datetime.date.today().isoformat()
summary.to_excel(
    f'/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN807/assays/BCR_sequencing/misc_files/data_processing/HVTN807_BCR_0428_upload_sample_summary_{today}.xlsx',
    index=True
)

manifest = pd.read_excel(
    '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN807/assays/BCR_sequencing/misc_files/manifests/HVTN 807 PBMC Shipped to FHCC as of 17Apr2026.xlsx'
)
manifest['guspec_core'] = manifest['Original Id'].str.rpartition("-")[0]
assert set(df_check.guspec_core).difference(manifest.guspec_core) == set()