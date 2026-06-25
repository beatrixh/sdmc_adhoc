## ---------------------------------------------------------------------------##
# Author: Sara Thiebaud
# Date: 24 June 2026
# Purpose: HVTN 300 blocking ELISA data from the Haynes lab; LDMS checks only
## ---------------------------------------------------------------------------##

import pandas as pd
import numpy as np
import sdmc_tools.access_ldms as access_ldms
import datetime
import os
import pdb
import datetime

INPUT_DATA_PATH = '/trials/vaccine/p300/s001/qdata/LabData/ELISA_blocking_pass-through/uploaded_by_lab/20260618-03/HVTN 300 Blocking Master File_v3.xlsx'


# LDMS CHECK ---------------------------------------------------------------------------------------------------- #
# read in data


ds = pd.read_excel(INPUT_DATA_PATH, skiprows=4)



COLUMN_NAMING = {
    'Subject ID':'ptid',
    'Visit':'visitno',
    'Original ID':'guspec',
    'Date Drawn':'drawdt'
}

ds = ds.rename(columns=COLUMN_NAMING)

#drop the non-record row (it's a comment)
ds.dropna(inplace=True, subset=['guspec'])

# pull ldms
ldms = access_ldms.pull_one_protocol('hvtn', 300)

# subset ldms
ldms = ldms.loc[ldms.guspec.isin(ds.guspec.dropna().unique())]

# check ptid and visit
m = ds.merge(ldms, on='guspec', how='left', suffixes=("_lab",""), indicator=True)
m.txtpid = m.txtpid.astype(float)
assert len(m.loc[(m.ptid!=m.txtpid) & (m.ptid.notna()),['ptid','txtpid']]) == 0

m.vidval = m.vidval.astype(float)
#strip the asterisks from visit numbers; the lab used these asterisks to denote reruns.
m['visitno'] = m.apply(lambda row: int(str(row['visitno']).strip('*')), axis=1)
# m['discrep_visitno'] = m.apply(lambda row: True if row['visitno']!=row['vidval'] else False, axis=1)
# m[m['discrep_visitno']==True].to_csv('LDMS_visitno_discrepancies_2026-06-24.csv', index=False)
assert len(m.loc[(m.visitno!=m.vidval) & (m.visitno.notna()),['visitno','vidval']]) == 0

#m['drawdt_ldms'] = m.apply(lambda row: datetime.date(int(row['drawdy']), int(row['drawdm']), int(row['drawdd'])) if not pd.isnull(row['drawdy']) else '', axis=1)
m['drawdt_d'] = m.apply(lambda row: row['drawdt'].day, axis=1)
m['drawdt_m'] = m.apply(lambda row: row['drawdt'].month, axis=1)
m['drawdt_y'] = m.apply(lambda row: row['drawdt'].year, axis=1)


assert len(m.loc[(m.drawdt_d!=m.drawdd) & (m.drawdt_d.notna()),['drawdt_d','drawdd']]) == 0
assert len(m.loc[(m.drawdt_m!=m.drawdm) & (m.drawdt_m.notna()),['drawdt_m','drawdm']]) == 0
assert len(m.loc[(m.drawdt_y!=m.drawdy) & (m.drawdt_y.notna()),['drawdt_y','drawdd']]) == 0