## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 2/23/2026
# Purpose: MAAR (process second upload with remaining guspecs and append to original)
## ---------------------------------------------------------------------------##
import pandas as pd
import numpy as np
import os
import datetime as datetime

import sdmc_tools.process as sdmc_tools
import sdmc_tools.access_ldms as access_ldms

# TODO: clean up notebook formatting

input_data_path = '/trials/vaccine/p206/s001/qdata/LabData/MAAR_pass-through/MAAR-Results-HVTN206-10Feb2026-Revised-20Feb2026.xlsx'
df = pd.read_excel(input_data_path)

df = df.rename(columns={'Global Spec IDs':'guspec'})

df = df.melt(
    id_vars=['PID','guspec'],
    value_vars=[
        'Abbott HIV-1/2 Ag/Ab',
        'BioRad HIV-1/2 Ag/Ab Combo',
        'BioRad HIV-1/2 Ab +O (3rd generation)',
        'Alere Determine HIV-1/2 Ag/Ab',
        'Diasorin Liaison HIV-1/2 Ag/Ab',
        'BioRad Geenius HIV-1/2 Ab',
        'INSTI HIV-1/2 Ab Rapid',
        'OraQuick HIV-1/2 Ab Rapid',
        'Abbott HIV-1 RNA',
    ],
    var_name='assay_subtype_lab',
    value_name='result'
)

df.loc[df.result.isna()]

result_detail = df.result.str.split(",", expand=True).rename(
    columns={0:'result_qualitative', 1:'result_quantitative'}
)

# old = pd.read_csv('/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN206/assays/MAAR/misc_files/data_processing/HVTN206_MAAR_processed_2026-01-22.txt', sep="\t")

df = pd.concat([
    df.drop(columns='result'),
    result_detail
], axis=1)

#move values in qual column into an 'result_detail' column for biorad geenius
df.loc[(df.assay_subtype_lab=='BioRad Geenius HIV-1/2 Ab') & (df.result_qualitative!='Non-reactive'), 'result_detail'] = df.loc[(df.assay_subtype_lab=='BioRad Geenius HIV-1/2 Ab') & (df.result_qualitative!='Non-reactive'), 'result_qualitative']

df

# amend biorad geenius qual column to contain "Ab Reactive" for all rows that aren't non-reactive
df.loc[(df.assay_subtype_lab=='BioRad Geenius HIV-1/2 Ab') & (df.result_qualitative.isin(['HIV-1 Ind. gp160', 'HIV Ind. gp140 and gp160'])), 'result_qualitative'] = 'Ab Reactive'



df.loc[(df.assay_subtype_lab=='BioRad Geenius HIV-1/2 Ab'),[
    'result_qualitative','result_quantitative','result_detail'
]].drop_duplicates()

# add units
df.loc[(df.assay_subtype_lab!='Abbott HIV-1 RNA') & (df.result_quantitative.notna()), 'result_units'] = 's/co'

# remove units from quant column
df.result_quantitative = df.result_quantitative.str.replace("s/co=","")

# add units cont'd
df.loc[(df.result_quantitative.notna()) & (df.result_quantitative.str.contains("cp/mL")), 'result_units'] = 'cp/mL'
df.result_quantitative = df.result_quantitative.str.replace("cp/mL","")
df.result_quantitative = df.result_quantitative.str.strip()

df.loc[(df.assay_subtype_lab=='Abbott HIV-1 RNA') & (df.result_qualitative!='Not Detected'), 'result_qualitative'] = 'Not Done'
df.loc[(df.assay_subtype_lab=='Abbott HIV-1 RNA') & (df.result_qualitative!='Not Detected'), 'result_quantitative'] = 'Not Done'

maar_metadata = pd.read_excel('/networks/vtn/lab/SDMC_labscience/assays/MAAR/UW/SDMC_materials/MAAR_info.xlsx')
maar_metadata = maar_metadata.iloc[:9]
maar_metadata = maar_metadata.rename(columns={'Name in incoming data':'assay_subtype_lab'})

metadata_usecols = [
    'assay_name',
    'assay_name_HAWS',
    'assay_subtype',
    'assay_precision',
    'instrument',
    'instrument_serial',
    'lab_software',
    'assay_subtype_lab',
    'LLOD',
]

maar_metadata.LLOD = maar_metadata.LLOD.str.replace("copies/mL","").str.strip()

df = df.merge(
    maar_metadata[metadata_usecols], on='assay_subtype_lab', how='outer'
)

df = df.drop(columns='assay_subtype_lab')
df = df.rename(columns={'instrument_serial':'instrument_serialno'})

ldms = access_ldms.pull_one_protocol('hvtn', 206)

md = {
    'network':'HVTN',
    'specrole':'Sample',
    'upload_lab_id':'UW',
    'assay_lab_name':'University of Washington Virology'
}

outputs = sdmc_tools.standard_processing(
    input_data=df,
    input_data_path=input_data_path,
    guspec_col='guspec',
    network='hvtn',
    metadata_dict=md,
    ldms=ldms,
)

assert (outputs.ptid.astype(int)!=outputs.pid.astype(int)).sum()==0
outputs.ptid = outputs.ptid.astype(int)

outputs = outputs.drop(
    columns=['pid']
)

reorder = [
    'network',
    'protocol',
    'guspec',
    'specrole',
    'upload_lab_id',
    'assay_lab_name',
    'assay_name_haws',
    'ptid',
    'visitno',
    'drawdt',
    'spectype',
    'spec_primary',
    'spec_additive',
    'spec_derivative',
    'assay_name',
    'assay_subtype',
    'assay_precision',
    'lab_software',
    'instrument',
    'instrument_serialno',
    'result_qualitative',
    'result_quantitative',
    'result_detail',
    'result_units',
    'llod',
    'sdmc_processing_datetime',
    'sdmc_data_receipt_datetime',
    'input_file_name',
]

set(reorder).symmetric_difference(outputs.columns)

outputs = outputs[reorder]

outputs = outputs.loc[~(
    (outputs.result_qualitative.isna()) & (outputs.result_quantitative.isna()) & (outputs.result_detail.isna())
)]

outputs.protocol= outputs.protocol.astype(int)
outputs.visitno= outputs.visitno.astype(float)
outputs.llod= outputs.llod.astype(float)

# tmp = outputs.copy()
# tmp['guspec_core'] = tmp.guspec.str.rpartition("-")[0]
# tmp = tmp.drop(columns='guspec')
# tmp = tmp.drop_duplicates()

# tmp.shape, prev_processed.shape

# tmp = tmp.loc[~(tmp.result_quantitative=="Not Done")]
# prev_tmp = prev_processed.loc[~(prev_processed.result_quantitative=="Not Done")]
# prev_tmp = prev_tmp.drop(columns='guspec_aliquots_possible')

# # tmp.shape, prev_tmp.shape

# set(tmp.columns).symmetric_difference(prev_tmp.columns)

# sort = tmp.columns.tolist()
# prev_tmp = prev_tmp[sort].sort_values(by=sort).reset_index(drop=True)
# tmp = tmp[sort].sort_values(by=sort).reset_index(drop=True)

# diffs = tmp.compare(prev_tmp)

# diffs['instrument_serialno'].drop_duplicates()

today = datetime.date.today().isoformat()
savedir = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN206/assays/MAAR/misc_files/data_processing/'
outputs.to_csv(savedir + f"HVTN206_MAAR_second_half_processed_{today}.txt", sep='\t', index=False)

prev_outputs = pd.read_csv(
    '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN206/assays/MAAR/misc_files/data_processing/HVTN206_MAAR_processed_2026-01-22.txt',
    sep='\t'
)

full_outputs = pd.concat([prev_outputs, outputs])

full_outputs[['result_qualitative','result_detail','result_units']].drop_duplicates()

full_outputs.to_csv(savedir + f"HVTN206_MAAR_full_dataset_processed_{today}.txt", sep='\t', index=False)

new_summary = outputs.loc[outputs.result_qualitative!='Not Done'].pivot(
    index=['ptid','visitno'],
    columns='assay_subtype',
    values='result_qualitative'
).fillna("Not Done")

new_summary.to_excel(
    savedir + f"HVTN206_MAAR_lab_feb_upload_summary_{today}.xlsx"
)

prev_full = pd.read_csv(
    '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN206/assays/MAAR/misc_files/data_processing/HVTN206_MAAR_full_dataset_processed_2026-02-13.txt',
    sep='\t')

prev_full[['result_qualitative','result_detail','result_units']].drop_duplicates()

full_outputs.shape == prev_full.shape

prev_full = prev_full[full_outputs.columns.tolist()]

compare = full_outputs.reset_index(drop=True).compare(prev_full.reset_index(drop=True))

compare.instrument_serialno.drop_duplicates()

compare.result_detail.drop_duplicates()

manifest = pd.read_csv(
    '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN206/assays/MAAR/misc_files/manifests/512-015-2025000266.txt',
    sep="\t"
)

set(outputs.guspec).difference(manifest.GLOBAL_ID)

set(outputs.guspec).difference(manifest.loc[manifest.PROTOCOL==206.].GLOBAL_ID)

set(manifest.loc[manifest.PROTOCOL==206.].GLOBAL_ID).difference(outputs.guspec)

50*9

outputs.loc[outputs.result_qualitative!='Not Done'].shape

set(outputs.guspec).symmetric_difference(manifest.loc[manifest.PROTOCOL==206.].GLOBAL_ID)


set(outputs.guspec).difference(lab_manifest['Global Spec ID'])
