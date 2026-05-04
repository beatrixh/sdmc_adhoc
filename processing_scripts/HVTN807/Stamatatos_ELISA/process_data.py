## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 05/01/2026
# Purpose:  run ad hoc processing + ldms checks
# NOTE -- this is quite drafty, will refine when i have metadata and confirmation of what some of those columns are
## ---------------------------------------------------------------------------##
import pandas as pd
import numpy as np
import sdmc_tools.access_ldms as access_ldms
import sdmc_tools.process as sdmc
import datetime
import os

input_data_path = "/trials/vaccine/p807/s001/qdata/LabData/ELISA_pass-through/HVTN807 ELISA data.xlsx"

df_all = pd.read_excel(input_data_path, sheet_name=None)

df = pd.DataFrame()
for k, tmp in df_all.items():
    if k!='ELISA Controls':
        tmp = tmp.dropna(how='all')
        tmp.columns = tmp.iloc[0]
        
        tmp = tmp.iloc[1:].melt(
            id_vars='Dilution',
            var_name='guspec',
            value_name='result'
        )
        tmp['k'] = k
        df = pd.concat([df, tmp])

# split data out into a few columns
df = pd.concat([
    df.k.str.replace(")","").str.split("(", expand=True).rename(columns={0:'group',1:'visitno'}),
    df.drop(columns='k')
], axis=1)

df.visitno = df.visitno.str.split(" ", expand=True)[0].astype(int)

df = pd.concat([
    df.group.str.strip().str.rpartition(" ", expand=True)[[0,2]].rename(columns={0:'group',2:'analyte_alternate'}),
    df.drop(columns='group')
], axis=1)

df = pd.concat([
    df.guspec.str.rpartition("\xa0", expand=True)[[0,2]].rename(columns={0:'guspec',2:'virus_dilution_maybe'}),
    df.drop(columns='guspec')
], axis=1)

df = pd.concat([
    df.guspec.str.rpartition("\xa0", expand=True)[[0,2]].rename(columns={0:'analyte',2:'guspec'}),
    df.drop(columns='guspec')
], axis=1)

df = df.loc[df.result.notna()]

df.virus_dilution_maybe = df.virus_dilution_maybe.str[1:-1]

# standard processing
ldms = access_ldms.pull_one_protocol('hvtn', 807)

md = {
    'network':'HVTN',
    'specrole':'Sample',
    'upload_lab_id':'I6',
    'assay_lab_name':'Stamatatos Lab',
    'assay_type':'ELISA',
    'assay_subtype':'',
    'assay_details':'',
    'assay_precision':'',
    'instrument':'',
    'lab_software':'',
    'lab_software_version':'',
    'result_units':'',
}

outputs = sdmc.standard_processing(
    input_data=df.drop(columns=['visitno']),
    input_data_path=input_data_path,
    guspec_col='guspec',
    network='HVTN',
    metadata_dict=md,
    ldms=ldms
)

today = datetime.date.today().isoformat()
savedir = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN807/assays/ELISA/misc_files/data_processing/'

outputs.to_csv(
    savedir+f"HVTN807_ELISA_processed_{today}.txt",
    sep="\t"
)

# pivot summary
summary = pd.pivot_table(
    df,
    index='txtpid',
    columns=['analyte'],
    aggfunc='count',
    values='guspec'
)

summary.to_excel(
    savedir+f'HVTN807_ELISA_ptid_analyte_sample_summary_{today}.xlsx'
)



# CHECK VISITNO / GROUP
df = df.merge(ldms[['guspec','vidval','txtpid']], on='guspec', how='left')
assert len(df.loc[df.visitno!=df.vidval]) == 0


# check group assignment against projected visits
pvisits = pd.read_excel(
    '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN807/assays/BCR_sequencing/misc_files/lab_temp_374998_2026-04-28_15-18-23.xlsx'
)

pvisits = pvisits[['PTID','Group']]
pvisits = pvisits.rename(columns = {'PTID':'txtpid','Group':'pvisits_group'})

df.txtpid = df.txtpid.astype(int)
pvisits.txtpid = pvisits.txtpid.astype(int)

pvisits.pvisits_group = pvisits.pvisits_group.str.rpartition(" ")[0]
df = df.merge(pvisits, on='txtpid', how='left')

# visual check -- these should be the same
df[['group','pvisits_group']].drop_duplicates()