import pandas as pd
import numpy as np
import sdmc_tools.access_ldms as access_ldms
import datetime
import os

# LDMS CHECK ---------------------------------------------------------------------------------------------------- #
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

# COMPLETENESS CHECK ---------------------------------------------------------------------------------------------------- #
pvisits = pd.read_excel(
    "/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN300/assays/lab_temp_331688_2026-06-23_10-36-18.xlsx"
)

pvisits = pvisits[['PTID','Note 11','Part']]
pvisits['in_data'] = pvisits.PTID.astype(int).isin(ldms.txtpid.astype(int))

# note, for Part B ptids, we only expect visit 11. per tracker entry:
    #EMPEM PT Report #2 (lab analysis): EMPEM will be performed by negative stain 
    #microscopy polyclonal epitope mapping (nsEMPEM) on serum (SST) samples from Part B V11 [M12.5, 2 wk. post 5th vacc.]

# Part B ptids NOT in data are all terminated
assert len(pvisits.loc[(pvisits.Part=="B") & (pvisits.in_data==False) & (pvisits['Note 11']!="Terminated")]) == 0

# all other part B ptids ARE in data
assert len(pvisits.loc[(pvisits.in_data==True) & (pvisits.Part!="B")]) == 0