import pandas as pd
import numpy as np
import sdmc_tools.access_ldms as access_ldms
import sdmc_tools.process as sdmc
import datetime
import os


input_data_path = "/trials/vaccine/p301/s001/qdata/LabData/BCR_sequencing_pass-through/uploaded_by_lab/McElrath/20260624-04/2026-06-24_HVTN301Boost_BCS2104-BCS2120_AIRR_filtered.tsv"
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
in_data_weird_in_pvisits = set(summary.reset_index().PTID.tolist()).difference(pvisits.loc[pvisits['Note 10']=='Complete'].PTID)

# ptids in pvisits but not the data
in_pvisits_not_data = set(pvisits.loc[pvisits['Note 10']!='Not Applicable'].PTID).difference(summary.reset_index().PTID)

# print to see
print(pvisits.loc[pvisits.PTID.isin(in_pvisits_not_data),['PTID','Visit 10','Note 10']])
print(pvisits.loc[pvisits.PTID.isin(in_data_weird_in_pvisits),['PTID','Visit 10','Note 10']])

# no sample for these, so dont expect them in the data
print(pvisits.loc[pvisits.PTID.isin(in_pvisits_not_data),['PTID','Visit 10','Note 10']])

# looks like a bunch that were out of window were included, and one that was NaN was included.
print(pvisits.loc[pvisits.PTID.isin(in_data_weird_in_pvisits),['PTID','Note 10','Note 14']])

# these two ptids shouldve had a 14
binary_summary = (summary > 0).astype(bool)
print(binary_summary[binary_summary[10] & ~binary_summary[14]])

questions = binary_summary[binary_summary[10] & ~binary_summary[14]].reset_index().PTID.tolist()

discrep = pvisits.loc[pvisits.PTID.isin(questions),['PTID','Note 10','Note 14']].rename(columns={'Note 10':'Visit 10', 'Note 14':'Visit 14'})
discrep = discrep.melt(id_vars='PTID', var_name='Visit', value_name='Projected Visits Status')

discrep.loc[discrep.Visit=='Visit 10', 'In Data'] = True
discrep.loc[discrep.Visit=='Visit 14', 'In Data'] = False

savedir="/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN301/assays/BCR_Sequencing/misc_files/data_processing/"
today = datetime.date.today().isoformat()

with pd.ExcelWriter(savedir+f'HVTN301_BCR_updated_sample_summary_{today}.xlsx', engine="openpyxl") as writer:
    binary_summary.to_excel(writer, sheet_name="summary", index=True)
    discrep.to_excel(writer, sheet_name="ptid_discrepancies", index=False)

# how does this compare to the may upload?
og = pd.read_csv(
    "/trials/vaccine/p301/s001/qdata/LabData/BCR_sequencing_pass-through/uploaded_by_lab/McElrath/20260518-02/2026-05-18_HVTN301_BCS2104-BCS2112_AIRR_filtered.tsv", 
    sep="\t"
)

# all guspecs from the old are in the new
assert set(og.Global_Spec_ID).difference(df.Global_Spec_ID) == set()

new_guspecs = set(df.Global_Spec_ID).difference(og.Global_Spec_ID)

# here are two guspecs in the new that werent in the old
print(ldms.loc[ldms.guspec_core.isin(set(df.Global_Spec_ID).difference(og.Global_Spec_ID))])

# unique row id
df['mid'] = df['cell_id'] + df['Global_Spec_ID']
og['mid'] = og['cell_id'] + og['Global_Spec_ID']

assert len(df) == df.mid.nunique()
assert len(og) == og.mid.nunique()

# things in old not in new
len(
    set(og.mid).difference(df.mid)
)

# most things in old are in new
len(
    set(og.mid).intersection(df.mid)
)

# cell_id/guspecs in the old but not the new. sara confirmed with nathifa over email, this was expected
missings = og.loc[og.mid.isin(set(og.mid).difference(df.mid)),['cell_id','Global_Spec_ID']]

missings.to_csv(
    '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN301/assays/BCR_Sequencing/misc_files/data_processing/in_may_upload_not_in_june_upload.csv',
    index=False
)


# look at the missing samples
problems = [777600479,821346268,726168211,764588462]
problems = pvisits.loc[pvisits.PTID.isin(problems),['PTID','Note 10','Note 14']]

problems.columns = ['PTID','10','14']

problems = problems.melt(id_vars='PTID', var_name='Visit', value_name='Status')

to_merge = df[['PTID','Visit']]
to_merge.Visit = to_merge.Visit.astype(str)
to_merge['in_data'] = True
problems = problems.merge(to_merge.drop_duplicates(), on=['PTID','Visit'], how='left')
problems=problems.fillna(False)

# it looks like the only ones that mightve been in the data but werent are:

# (764588462/14): sample 10 doesnt exist, shouldnt be included 
# (777600479/14) sample 10 was "NA", maybe shouldnt be included? 
# (821346268/14): sample 10 was out of window, maybe shouldnt be included?


# for all the others the sample doesnt exist