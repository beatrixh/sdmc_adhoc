## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 04/30/2026
# Purpose:  Verify that the wide and the long BCR data looks identical
# Note:		We only added columns (didnt touch any of ju's outputs at all), so can do this on our outputs.
## ---------------------------------------------------------------------------##

import pandas as pd
import numpy as np
import os

# read in data
savedir = "/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN807/assays/BCR_sequencing/misc_files/data_processing/"
wide = pd.read_csv(
    savedir+"HVTN807_BCR_sequencing_regenerated_AIRR_wide_processed_2026-04-29.txt",
    sep="\t"
)

long = pd.read_csv(
    savedir+"HVTN807_BCR_sequencing_regenerated_AIRR_long_processed_2026-04-29.txt",
    sep="\t"
)


# drop file-specific processing stamps that we added
wide = wide.drop(columns=['input_file_name','sdmc_processing_datetime','sdmc_data_receipt_datetime'])
long = long.drop(columns=['input_file_name','sdmc_processing_datetime','sdmc_data_receipt_datetime'])

# these are the columns shared between the datasets
shared = list(
    set(wide.columns).intersection(long.columns)
)

# can we convert the wide to the long
wide_melt = wide.melt(
    id_vars=shared,
    value_vars=list(set(wide.columns).difference(shared))
)

wide_melt = pd.concat([
    wide_melt.drop(columns='variable'),
    wide_melt.variable.str.rpartition("_", expand=True)[[0,2]].rename(columns={0:'variable', 2:'chain'})
], axis=1)

wide_to_long = wide_melt.pivot(
    index=shared + ['chain'],
    columns='variable',
    values='value'
)

# make sure these have the same column and row ordering
wide_to_long= wide_to_long[list(long.columns)]
wide_to_long = wide_to_long.sort_values(by=['guspec_core','cell_id','chain']).reset_index(drop=True)
long = long.sort_values(by=['guspec_core','cell_id','chain']).reset_index(drop=True)

# this is the one difference -- looks like this accidentally has some float-valued 1s instead of string-valued-int-shaped?
long.loc[long.v2apex_candidate_type==1,'v2apex_candidate_type'] = '1'

# check that we successfully converted the wide to be identical to the long
compare = wide_to_long.compare(long)
assert len(compare) == 0

