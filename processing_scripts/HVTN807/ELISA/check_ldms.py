import pandas as pd
import numpy as np
import sdmc_tools.access_ldms as access_ldms
import sdmc_tools.process as sdmc
import datetime
import os

input_data_path = "/trials/vaccine/p807/s001/qdata/LabData/ELISA_pass-through/HVTN807 ELISA data.xlsx"
df_all = pd.read_excel(input_data_path, sheet_name=None)

gmap = pd.DataFrame()
for k, v in df_all.items():
    guspecs = v.iloc[1,1:].dropna().tolist()
    s = pd.DataFrame(data={'k':[k]*len(guspecs), 'guspecs':guspecs})
    gmap = pd.concat([gmap, s])

# dont need to check these
df = gmap.loc[gmap.k!="ELISA Controls"]

# some formatting -- need to pull out guspecs, visitnos, and groups
df = pd.concat([
    df.k.str.replace(")","").str.split("(", expand=True).rename(columns={0:'group',1:'visitno'}),
    df[['guspecs']]
], axis=1)

df.visitno = df.visitno.str.split(" ", expand=True)[0].astype(int)

df = pd.concat([
    df.group.str.strip().str.rpartition(" ", expand=True)[[0,2]].rename(columns={0:'group',2:'something'}),
    df.drop(columns='group')
], axis=1)

df = pd.concat([
    df.guspecs.str.rpartition("\xa0", expand=True)[[0,2]].rename(columns={0:'guspec',2:'something2'}),
    df.drop(columns='guspecs')
], axis=1)

df = pd.concat([
    df.guspec.str.rpartition("\xa0", expand=True)[[0,2]].rename(columns={0:'something3',2:'guspec'}),
    df.drop(columns='guspec')
], axis=1)


# check visitno against ldms
ldms = access_ldms.pull_one_protocol('hvtn', 807)
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