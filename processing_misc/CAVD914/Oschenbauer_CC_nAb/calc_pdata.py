## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 2025-10-31
# Purpose: Process CC nAb data from Ochsenbauer
# 			- outputs: mab_name, isolate, pubid, ic50, ic80
#			- use curve estimates from nab tool
#			- when multiple replicates, take geometric mean
#			- for now, drop anything problematic or ambiguous
## ---------------------------------------------------------------------------##

import pandas as pd
import numpy as np
import os
import scipy.stats as stats
import datetime

# read in data
df = pd.read_csv('/networks/cavd/Studies/cvd914/qdata/LabData/NAB/CAVD914_NAB07_U6_20250924.txt', sep='\t')

# reformat columns
df.columns = [i.lower() for i in df.columns]
df = df.rename(columns={
    'isolate':'tetO_HIV_Env_IMC_name',
    'guspec':'mab_name',
})
# pull isolate out
df['isolate'] = df.tetO_HIV_Env_IMC_name.str.split(".", expand=True)[1]

# pull pubID out
pubids = df.isolate.str.split("_", expand=True).iloc[:,:2]
df['pubid'] = pubids[0]+"-"+pubids[1]

# check one isolate per pubID
assert set(df.groupby('pubid').isolate.nunique())=={1}

# subset to IC50/80s + ids
ics = df.loc[df.value_label.isin(['TITER_50_CURVE','TITER_80_CURVE']),['mab_name','isolate','pubid','value_label','value']]
ics['id_var'] = ics.mab_name+"|"+ics.isolate+"|"+ics.value_label

# filter out anything that was below or above nAb tool cutoff
ics['above_or_below_cutoff'] = False
ics.loc[ics.value.str.contains(">"), 'above_or_below_cutoff'] = True
ics.loc[ics.value.str.contains("<"), 'above_or_below_cutoff'] = True
filter_ = ics.groupby('id_var').above_or_below_cutoff.sum()
filter_ = filter_.loc[filter_==0].reset_index().id_var.tolist()
unambig = ics.loc[ics.id_var.isin(filter_)]

# should only contain unambiguous mab-isolate pairs, now
assert unambig.above_or_below_cutoff.sum() == 0

# calculate geometric mean
unambig.value = unambig.value.astype(float)
unambig = unambig.groupby(['mab_name','isolate','pubid','value_label'])['value'].apply(stats.gmean).reset_index()

# cast to wide along IC50 v 80
unambig.value_label = unambig.value_label.map({
    'TITER_50_CURVE':'ic50',
    'TITER_80_CURVE':'ic80',
})
outputs = unambig.pivot(
    index=['mab_name','isolate','pubid'],
    columns='value_label',
    values='value'
).reset_index()


# is there anything weird
(outputs.ic50 >= outputs.ic80).sum()

# yes there is anything weird, get rid of it for now
outputs_sensible = outputs.loc[outputs.ic50 < outputs.ic80]

# save to .txt
savedir = '/networks/vtn/lab/SDMC_labscience/studies/VISC/Ochsenbauer_914/assays/C-C_nAb/misc_files/example_data/'
today = datetime.date.today().isoformat()
outputs_sensible.to_csv(savedir + f"CAVD914_Ochsenbauer_CC_nAb_subset_example_{today}.txt", sep='\t', index=False)