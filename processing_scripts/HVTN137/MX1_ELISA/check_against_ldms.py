import pandas as pd
import numpy as np
import sdmc_tools.access_ldms as access_ldms
import sdmc_tools.process as sdmc
import os

datadir = '/trials/vaccine/p137/s001/qdata/LabData/MX1_ELISA_pass-through/'
fnames = [i for i in os.listdir(datadir) if 'metadata' in i and '$' not in i]

data = pd.concat([pd.read_excel(datadir + f) for f in fnames])

# subset to rows in LDMS
for i in data.loc[data.guspec.isna()].specrole.unique():
    assert i in ['QC Control', 'Standard']

data = data.loc[data.guspec.notna()]

# subset to the cols we want to check
checkcols = [
    'network',
    'protocol',
    'guspec',
    'ptid',
    'drawdt',
    'visitno',
    'spectype',
]
data = data[checkcols]

data.columns = [i + "_og" for i in data.columns]

# 0373-01XTRF00-009 has a leading space, need to clean for LDMS merge
data.loc[data.guspec_og.str.contains(" ")].guspec_og.tolist()
data.guspec_og = data.guspec_og.str.strip()

# merge on lDMS
ldms = access_ldms.pull_one_protocol('hvtn', 137)
outputs = sdmc.add_ldms(
    ldms,
    data,
    guspec_col='guspec_og',
)

# LDMS CHECKS
assert (outputs.visitno.astype(float)!=outputs.visitno_og.astype(float)).sum() == 0
assert (outputs.ptid.astype(int)!=outputs.ptid_og.astype(int)).sum() == 0
assert (outputs.protocol.astype(int)!=outputs.protocol_og.astype(int)).sum() == 0

# outputs.drawdt = pd.to_datetime(outputs.drawdt).dt.date
outputs.drawdt_og = pd.to_datetime(outputs.drawdt_og).dt.date
assert (outputs.drawdt.astype(str)!=outputs.drawdt_og.astype(str)).sum() == 0

# this was saved as chilled serum, whereas LDMS has serum. this is ok.
print(outputs.loc[outputs.spectype!=outputs.spectype_og,['spectype','spectype_og']].drop_duplicates())