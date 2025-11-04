## ------------------------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 2025-11-03
# Purpose: Process CC nAb data from Ochsenbauer
#           - outputs: mab_name, isolate, pubid, ic50, ic80
#           - use curve estimates from nab tool
#           - discard rows for which nab tool didnt fit data
#           - when multiple replicates, take geometric mean
#           - when at least one replicate has a censored value, drop other replicates
## ------------------------------------------------------------------------------------------##

import pandas as pd
import numpy as np
import os
import scipy.stats as stats
import datetime

## pull out the rows that the nAb tool didn't fit (we want to exclude these) ----------------##
atlas = pd.read_csv(
    '/networks/vtn/lab/SDMC_labscience/studies/VISC/Ochsenbauer_914/assays/C-C_nAb/misc_files/atlas_export/NAb_CAVD_OCH_2025-10-31_15-35-11.tsv',
    sep='\t'
)

atlas.ISOLATE = atlas.ISOLATE.str.split(".", expand=True)[1].str[1:]
atlas = atlas.rename(columns={'SPECID':'mab_name'})
atlas.columns = [i.lower().replace(" ","_") for i in atlas.columns]
atlas.percent_neutralization = atlas.percent_neutralization.str[:-1].astype(float)

atlas = atlas.rename(columns={
    'mean':'sample_mean',
    'std_dev':'sample_std_dev',
})

# pull rows where sample RLU is NaN
no_fits = (atlas.loc[atlas.sample_mean.isna()].mab_name + "|" + atlas.loc[atlas.sample_mean.isna()].isolate).unique().tolist()

# read in data with IC50/80 estimates -------------------------------------------------------##
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

# confirm one isolate per pubID
assert set(df.groupby('pubid').isolate.nunique())=={1}

df['mab_isolate_titer'] = df['mab_name']+"|"+df['isolate']+"|"+df['value_label']
df['mab_isolate_pubid_titer'] = df['mab_name']+"|"+df['isolate']+"|"+df['pubid']+"|"+df['value_label']

ics = df.loc[df.value_label.isin(['TITER_50_CURVE','TITER_80_CURVE']),['mab_isolate_pubid_titer', 'mab_isolate_titer','mab_name','isolate','pubid','value_label','val']]

# mark rows corresponding to a value above/below the cutoff threshold
ics['above_cutoff'] = False
ics.loc[ics.val.str.contains(">"), 'above_cutoff'] = True

ics['below_cutoff'] = False
ics.loc[ics.val.str.contains("<"), 'below_cutoff'] = True

# look at range of values
# check = ics.val.copy()
# check.loc[ics.above_cutoff] = 100
# check.loc[ics.below_cutoff] = -1
# check = check.astype(float)
# check = list(np.sort(check.astype(float).unique()))
# check[:5] + check[-5:]

ic_zeros = (ics.loc[(ics.val=='0.0')].mab_name + "|" + ics.loc[(ics.val=='0.0')].isolate).unique().tolist()

# confirm that the IC50/80 estimates that are exactly zero come from rows for which nab tool didnt fit a curve / nan-data
assert len(set(ic_zeros).symmetric_difference(no_fits)) == 0

# subset to values from actual estimates
ics = ics.loc[ics.val!='0.0']

# most simple version -- assume marginals are on the boundary -------------------------------##
def calc_gmean_with_boundary_exceptions(x):
    above = False
    below = False
    if len(x) == 1:
        return str(x.values[0])
    else:
        for xi in x:
            if ">" in xi:
                above = True
            elif "<" in xi:
                below = True
        if above and below:
            raise Exception("Surprising result")
        elif above or below:
            # return all unique values concatenated with ','s
            return ','.join(set(x))
        else:
            return str(stats.gmean(x.astype(float)))

simplest = ics.copy()
simplest = simplest.groupby(['mab_isolate_pubid_titer','mab_name','isolate','pubid','value_label']).val.apply(calc_gmean_with_boundary_exceptions).reset_index()

simplest.value_label = simplest.value_label.map({
    'TITER_50_CURVE':'ic50',
    'TITER_80_CURVE':'ic80',
})

simplest = simplest.pivot(
    index=['mab_name','isolate','pubid'],
    columns='value_label',
    values='val'
).reset_index()

# mark estimates for which some but not all replicates were above/below cutoff
simplest['uncensored_values_discarded'] = False
simplest.loc[simplest.ic50.str.contains(","),'uncensored_values_discarded'] = 'IC50'
simplest.loc[simplest.ic80.str.contains(","),'uncensored_values_discarded'] = 'IC80'
simplest.loc[(simplest.ic50.str.contains(",")) & (simplest.ic80.str.contains(",")),'uncensored_values_discarded'] = 'IC50 and IC80'

# push all such estimates to cutoff
simplest.loc[simplest.ic50.str.contains(",") & (simplest.ic50.str.contains(">25")), 'ic50'] = '>25'
simplest.loc[simplest.ic50.str.contains(",") & (simplest.ic50.str.contains("<")), 'ic50'] = "<0.011431184270690446"

simplest.loc[simplest.ic80.str.contains(",") & (simplest.ic80.str.contains(">25")), 'ic80'] = '>25'
simplest.loc[simplest.ic80.str.contains(",") & (simplest.ic80.str.contains("<")), 'ic80'] = "<0.011431184270690446"

# save to .txt ------------------------------------------------------------------------------##
savedir = '/networks/vtn/lab/SDMC_labscience/studies/VISC/Ochsenbauer_914/assays/C-C_nAb/misc_files/data_processing/'
today = datetime.date.today().isoformat()
simplest.to_csv(savedir + f"CAVD914_Ochsenbauer_CC_nAb_IC50_80s_{today}.txt", sep='\t', index=False)

# additional checks--------------------------------------------------------------------------##
checks = simplest.copy()

checks['ic50'] = checks.ic50.str.replace(">","").str.replace("<","").astype(float)
checks['ic80'] = checks.ic80.str.replace(">","").str.replace("<","").astype(float)

checks.loc[checks.ic50 >= checks.ic80,['ic50','ic80']].drop_duplicates()

assert (checks.ic50 == 0.).sum() == 0
assert (checks.ic80 == 0.).sum() == 0
assert checks.ic50.isna().sum() == 0
assert checks.ic80.isna().sum() == 0

checks.loc[checks.uncensored_values_discarded=="IC50"].ic50.unique()
checks.loc[checks.uncensored_values_discarded=="IC80"].ic80.unique()
checks.loc[checks.uncensored_values_discarded=="IC50 and IC80",['ic50','ic80']].drop_duplicates()

# this one a little weird, seems like the censored values were probably an outlier for one side or the other, at least
checks.iloc[507]


