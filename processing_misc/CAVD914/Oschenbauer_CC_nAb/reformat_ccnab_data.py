## ------------------------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 2025-11-14
# Purpose: Reformat IC50/80 estimates from R script
## ------------------------------------------------------------------------------------------##

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datetime
import os

titer_estimates = pd.read_csv("/networks/vtn/lab/SDMC_labscience/studies/VISC/Ochsenbauer_914/assays/C-C_nAb/misc_files/interim_data/tobit_titers_2025-11-17.txt", sep='\t')

titer_estimates['titer_numeric'] = titer_estimates.titer.str.replace(">","").str.replace("<","").astype(float)
titer_estimates['titer_recensored'] = titer_estimates.titer.copy()
titer_estimates.loc[titer_estimates.titer_numeric > 25, 'titer_recensored'] = ">25"
titer_estimates.loc[titer_estimates.titer_numeric < 0.011431184270690446, 'titer_recensored'] = '<0.011431184270690446'

# ensure only one estimate per isolate/mab/titer
assert len(titer_estimates.drop_duplicates().groupby(['isolate','mab_name','pct_neut']).titer.count().value_counts()) == 1

titer_estimates = titer_estimates.sort_values(by=['isolate','mab_name','pct_neut']).drop_duplicates()
titer_estimates.pct_neut = titer_estimates.pct_neut.map({'TITER_50_CURVE':'IC50','TITER_80_CURVE':'IC80'})

titer_estimates['n_replicates'] = titer_estimates.n_numeric + titer_estimates.n_above + titer_estimates.n_below
titer_estimates['n_outside_cutoff'] = titer_estimates.n_above + titer_estimates.n_below

# make sure same number of replicates for ic50 and ic80
tmp = titer_estimates.pivot(
    index=['isolate','mab_name'],
    columns='pct_neut',
    values='n_replicates'
)
assert (tmp.IC50 != tmp.IC80).sum() == 0

recensored_wide = titer_estimates.pivot(
    index=['isolate','mab_name'],
    columns='pct_neut',
    values='titer_recensored'
)

precensored_wide = titer_estimates.pivot(
    index=['isolate','mab_name'],
    columns='pct_neut',
    values='titer'
)

outside_cutoff = titer_estimates.pivot(
    index=['isolate','mab_name'],
    columns='pct_neut',
    values='n_outside_cutoff'
).rename(columns={'IC50':'n_outside_cutoff_IC50','IC80':'n_outside_cutoff_IC80'})

estimates_wide = recensored_wide.merge(precensored_wide, left_index=True, right_index=True, suffixes=('_recensored','_tobit_fit'))
estimates_wide = estimates_wide.merge(outside_cutoff, left_index=True, right_index=True)

n_reps = titer_estimates.set_index(['isolate','mab_name'])[['n_replicates']].reset_index().drop_duplicates().set_index(['isolate','mab_name'])

estimates_wide = estimates_wide.merge(
    n_reps, 
    left_index=True, right_index=True
)

estimates_wide['n_replicates_ic50_numeric'] = estimates_wide.n_replicates - estimates_wide.n_outside_cutoff_IC50
estimates_wide['n_replicates_ic80_numeric'] = estimates_wide.n_replicates - estimates_wide.n_outside_cutoff_IC80

estimates_wide.columns = [i.lower() for i in estimates_wide.columns]

rename = {
    'ic50_recensored': 'ic50_recensored',
    'ic80_recensored': 'ic80_recensored',
    'ic50_tobit_fit': 'ic50',
    'ic80_tobit_fit': 'ic80',
    'n_outside_cutoff_ic50': 'n_outside_cutoff_ic50',
    'n_outside_cutoff_ic80': 'n_outside_cutoff_ic80',
    'n_replicates': 'n_replicates',
    'n_replicates_ic50_numeric': 'n_replicates_ic50_numeric',
    'n_replicates_ic80_numeric': 'n_replicates_ic80_numeric'
}
estimates_wide = estimates_wide.rename(columns=rename)

final = estimates_wide[['n_replicates','ic50','ic50_recensored','n_replicates_ic50_numeric','ic80','ic80_recensored','n_replicates_ic80_numeric']].reset_index()

# read in data with IC50/80 estimates to get pubid back
df = pd.read_csv('/networks/cavd/Studies/cvd914/qdata/LabData/NAB/CAVD914_NAB07_U6_20250924.txt', sep='\t')
df.columns = [i.lower() for i in df.columns]
df = df.rename(columns={
    'isolate':'tetO_HIV_Env_IMC_name',
    'guspec':'mab_name',
    'value':'val',
})
df['isolate'] = df.tetO_HIV_Env_IMC_name.str.split(".", expand=True)[1].str[1:]

pubids = df.isolate.str.split("_", expand=True).iloc[:,:2]
df['pubid'] = pubids[0]+"-"+pubids[1]

# these should uniquely define one another
assert df[['pubid','isolate']].drop_duplicates().shape[0] == df.pubid.nunique()
assert df.pubid.nunique() == df.isolate.nunique()

df[['isolate','pubid']].drop_duplicates()

final = final.merge(df[['pubid','isolate']].drop_duplicates(), how='left', on='isolate')

final = final[['isolate','mab_name','pubid','n_replicates','ic50','ic50_recensored','n_replicates_ic50_numeric','ic80','ic80_recensored','n_replicates_ic80_numeric']]

today = datetime.date.today().isoformat()
cc_dir = "/networks/vtn/lab/SDMC_labscience/studies/VISC/Ochsenbauer_914/assays/C-C_nAb/misc_files/"
final.to_csv(cc_dir + f"data_processing/CAVD914_CC_nAb_IC50_80_{today}.txt", sep="\t", index=False)

ic50_numeric = final.ic50.str.replace(">","").str.replace("<","").astype(float)
ic80_numeric = final.ic80.str.replace(">","").str.replace("<","").astype(float)

# visual check on outliers
final.loc[(ic50_numeric > ic80_numeric)]