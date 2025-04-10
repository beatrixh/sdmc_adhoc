## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 04/09/2025
# Purpose:  output
#           1. table comparing new vs old resutls as submitted
#           2. scatterplot comparing new vs old results as submitted
## ---------------------------------------------------------------------------##
import pandas as pd
import numpy as np
import datetime
import matplotlib.pyplot as plt

import os

# read in data
datadir = "/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN302/assays/AE_assays/Anti-PEG_total/misc_files/data_processing/"

new = pd.read_csv(datadir + "HVTN302_Anti-PEG_total_CRL_interpolated_titers_processed_2025-04-09.txt", sep="\t")
old = pd.read_csv(datadir + "archive_2025-04-04/HVTN302_Anti-PEG_total_CRL_processed_2024-07-02.txt", sep="\t")

# table comparing old vs new results as submitted ----------------------------##

compare = new[['guspec','ptid','visitno','result_as_submitted']].merge(
    old[['guspec','result_as_submitted']],
    on='guspec',
    how='outer',
    suffixes=(
    '_new','_old')
)

compare = compare.rename(columns={'result_as_submitted_new':'interpolated_result_2025_04_09',
                        'result_as_submitted_old':'endpoint_titer_original_2024_07_02'})
                        compare.to_csv("/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN302/assays/AE_assays/Anti-PEG_total/misc_files/HVTN302_Anti-PEG_total_original_versus_rerun_interpolated_titers.txt", sep="\t")

# scatterplot ----------------------------------------------------------------##

# merge together
both = new[['guspec','result_as_submitted']].merge(
    old[['guspec','result_as_submitted']],
    on='guspec',
    how='outer',
    suffixes=(
    '_new','_old')
)

# map results columns to numeric values
numeric_map = {
    'Negative Immunodepletion': np.nan,
    'Negative Screen': np.nan,
    'Negative Titer (<50.000)': np.nan,
    'Positive Titer (>6400.000)': 6400.
}
for version in ['old','new']:
    numeric_label = f"result_numeric_{version}"
    submitted_label = f"result_as_submitted_{version}"
    both[numeric_label] = both[submitted_label]
    both.loc[both[numeric_label].isin(numeric_map.keys()), numeric_label] = both.loc[both[numeric_label].isin(numeric_map.keys()), numeric_label].map(numeric_map)
    both[numeric_label] = both[numeric_label].astype(float)



# separate out into rerun vs not-rerun results
not_rerun = both.loc[both.result_as_submitted_old!="Positive Titer (>6400.000)"]
rerun = both.loc[both.result_as_submitted_old=="Positive Titer (>6400.000)"]

# axis ticks
ticks = [50*2**i for i in range(10)]
yticks = ticks.copy()
xticks = ticks.copy()
xticks[-3:] = [">6400","",""]

# annotation
extra_text = 'Results with "Negative Titer\n(<50)" are not displayed\nhere. None of those\nresults changed'

# build plot
fig, ax = plt.subplots(figsize=(5,5))
lims = [min(ticks), max(ticks)]

# x=y line
ax.plot(lims, lims, '--', alpha=0.75, zorder=0, color='black', linewidth=1)

# titers
for titer in ticks:
    ax.axhline(titer, linewidth=0.3, color='gray')
    ax.axvline(titer, linewidth=0.3, color='gray')

# scatter data
ax.scatter(not_rerun.result_numeric_old, not_rerun.result_numeric_new, s=15, alpha=0.2, label='Not Rerun')
ax.scatter(rerun.result_numeric_old, rerun.result_numeric_new, s=15, alpha=0.2, label='Rerun')

# add notes
ax.annotate(extra_text, (1000,50))

# axes
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xticks(ticks, xticks, rotation=45)
ax.set_yticks(ticks, yticks)
ax.set_aspect('equal')
ax.set_xlabel("Original Endpoint Titers")
ax.set_ylabel("Interpolated Titers")
ax.legend()
plt.suptitle("Anti-PEG Result Updates")

savedir = "/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN302/assays/AE_assays/Anti-PEG_total/misc_files/"
plt.savefig(savedir + f"Anti-PEG_result_updates_scatterplot.png", dpi=320, format='png', transparent=False, bbox_inches='tight', pad_inches=0.3)
