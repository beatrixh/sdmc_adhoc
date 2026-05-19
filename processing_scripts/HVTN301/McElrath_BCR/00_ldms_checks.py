import pandas as pd
import numpy as np
import sdmc_tools.access_ldms as access_ldms
import sdmc_tools.process as sdmc
import datetime
import os


input_data_path = '/trials/vaccine/p301/s001/qdata/LabData/BCR_sequencing_pass-through/uploaded by lab/McElrath/2026-05-18_HVTN301_BCS2104-BCS2112_AIRR_filtered.tsv'
df = pd.read_csv(input_data_path, sep="\t")

# pivot summary
summary= df.groupby(['PTID','Visit']).count().max(axis=1).reset_index()
summary = summary.pivot(
    index='PTID',
    columns='Visit',
    values=0
).fillna(0)

# ldms checks
ldms = access_ldms.pull_one_protocol('hvtn', 301)
ldms['guspec_core'] = ldms.guspec.str.rpartition("-")[0]
ldms = ldms.drop(columns='guspec').drop_duplicates()

df_check = df[['Global_Spec_ID','Visit','PTID']].drop_duplicates().rename(columns={'PTID':'ptid_lab','Global_Spec_ID':'guspec_core','Visit':'visitno_lab'})
df_check.guspec_core = df_check.guspec_core.str.replace("\t","")

df_check = df_check.merge(ldms, on='guspec_core', how='left')
df_check = df_check.rename(columns={'txtpid':'ptid_ldms', 'vidval':'visitno_ldms'})

assert len(df_check.loc[(df_check.visitno_lab.astype(float) != df_check.visitno_ldms.astype(float)),['guspec_core','visitno_lab','visitno_ldms']]) == 0
assert len(df_check.loc[(df_check.ptid_lab.astype(float)!=df_check.ptid_ldms.astype(float)),['guspec_core','ptid_lab','ptid_ldms']]) == 0

# check completeness against projected visits
# we expect visit 10 and 14 for every ptid that had a visit 10
pvisits = pd.read_excel("/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN301/assays/lab_temp_44805_2026-05-19_11-27-35.xlsx")

# ptids in the data but not pvisits
in_data_not_pvisits = set(summary.reset_index().PTID.tolist()).difference(pvisits.loc[pvisits['Note 10']!='Not Applicable'].PTID)
in_data_not_pvisits

# ptids in pvisits but not the data
in_pvisits_not_data = set(pvisits.loc[pvisits['Note 10']!='Not Applicable'].PTID).difference(summary.reset_index().PTID)
in_pvisits_not_data

# print to see
print(pvisits.loc[pvisits.PTID.isin(in_pvisits_not_data),['PTID','Visit 10','Note 10']])

print(pvisits.loc[pvisits.PTID.isin(in_data_not_pvisits),['PTID','Visit 10','Note 10']])

ptid_discrepancies = pd.concat([
    pd.DataFrame({
        'ptid':list(in_data_not_pvisits),
        'note':["in data; 'Not Applicable' in projected visits"]
        }),
    pd.DataFrame({
        'ptid':[777340095],
        'note':["in projected visits, not in data"]
        }),
    ])

# a lot of discrepancies -- save to .xlsx
savedir="/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN301/assays/BCR_Sequencing/misc_files/data_processing/"
today = datetime.date.today().isoformat()

binary_summary = (summary > 0).astype(bool)
with pd.ExcelWriter(savedir+f'HVTN301_BCR_updated_sample_summary_{today}.xlsx', engine="openpyxl") as writer:
    binary_summary.to_excel(writer, sheet_name="summary", index=True)
    ptid_discrepancies.to_excel(writer, sheet_name="ptid_discrepancies", index=False)