import pandas as pd
import numpy as np
import sdmc_tools.access_ldms as access_ldms
import datetime
import os

# read in data
input_data_path = '/trials/vaccine/p300/s001/qdata/LabData/EMPEM_passthrough/HVTN300B_v11 NSEM Data 04May26.xlsx'
ptcls = pd.read_excel(input_data_path, sheet_name='Fab-bound_Ptcls')
igg = pd.read_excel(input_data_path, sheet_name='IgG Fab Recovery', skiprows=2)

ptcl_rename = {
    'NSEM ID':'lab_id1',
    'BSI ID':'lab_id2',
    'Subject ID':'ptid1',
    'Serum Original ID':'guspec',
    'PID':'ptid2',
}
ptcls = ptcls.rename(columns=ptcl_rename)

igg_rename = {
    'NSEM ID':'lab_id1',
    'BSI ID':'lab_id2',
    'Subject ID':'ptid1',
    'Serum Original ID':'guspec',
}
igg = igg.rename(columns=igg_rename)

# pull ldms
ldms = access_ldms.pull_one_protocol('hvtn', 300)


# verify two sheets have same guspecs
assert set(ptcls.guspec).difference(igg.guspec).difference(['Average', 'Maximum', 'Minimum', 'Standard Deviation']) == set()
assert set(igg.guspec.dropna()).difference(ptcls.guspec) == set()

# subset ldms
ldms = ldms.loc[ldms.guspec.isin(igg.guspec.dropna().unique())]

# check igg
igg = igg.merge(ldms, on='guspec', how='left', suffixes=("_lab",""))
igg.txtpid = igg.txtpid.astype(float)
assert len(igg.loc[(igg.ptid1!=igg.txtpid) & (igg.ptid1.notna()),['ptid1','txtpid']]) == 0

# check ptcls
ptcls = ptcls.merge(ldms, on='guspec', how='left', suffixes=("_lab",""))
ptcls = ptcls.loc[ptcls.guspec.notna()]
assert len(ptcls.loc[(ptcls.ptid1!=ptcls.txtpid.astype(float)) & (ptcls.ptid1.notna())]) == 0
assert len(ptcls.loc[(ptcls.ptid2!=ptcls.txtpid.astype(float)) & (ptcls.ptid2.notna())]) == 0

# all of them are visit 11
ldms.vidval.unique()