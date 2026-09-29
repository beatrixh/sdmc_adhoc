## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 2026-09-26
# Purpose: MAAR (process third upload with remaining guspecs and append to original)
## ---------------------------------------------------------------------------##
import pandas as pd
import numpy as np
import os
import datetime as datetime

import sdmc_tools.process as sdmc_tools
import sdmc_tools.access_ldms as access_ldms

input_data_path = '/trials/vaccine/p206/s001/qdata/LabData/MAAR_pass-through/MAAR-Results-HVTN206-24Sep2026.xlsx'
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

df.loc[df.result=="HIV-1 Pos. gp 41 and gp160, HIV-2 Ind. gp140", 'result'] = "HIV-1 Pos. gp 41 and gp160; HIV-2 Ind. gp140"

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

# amend biorad geenius qual column to contain "Ab Reactive" for all rows that aren't non-reactive
df.loc[(df.assay_subtype_lab=='BioRad Geenius HIV-1/2 Ab') & (df.result_qualitative.isin([
    'HIV-1 Ind. gp160',
    'HIV-2 Ind. gp140',
    'HIV Ind. gp140 and gp160',
    'HIV-1 Pos. gp 41 and gp160; HIV-2 Ind. gp140'
    
])), 'result_qualitative'] = 'Ab Reactive'

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

assert set(df.assay_subtype_lab.unique()).symmetric_difference(maar_metadata.assay_subtype_lab.unique()) == set()
df_w_metadata = df.merge(
    maar_metadata[metadata_usecols], on='assay_subtype_lab', how='outer'
)

df_w_metadata = df_w_metadata.drop(columns='assay_subtype_lab')
df_w_metadata = df_w_metadata.rename(columns={'instrument_serial':'instrument_serialno'})

ldms = access_ldms.pull_one_protocol('hvtn', 206)

md = {
    'network':'HVTN',
    'specrole':'Sample',
    'upload_lab_id':'UW',
    'assay_lab_name':'University of Washington Virology'
}

outputs = sdmc_tools.standard_processing(
    input_data=df_w_metadata,
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

assert set(reorder).symmetric_difference(outputs.columns) == set()

outputs = outputs[reorder]

outputs = outputs.loc[~(
    (outputs.result_qualitative.isna()) & (outputs.result_quantitative.isna()) & (outputs.result_detail.isna())
)]

outputs.protocol= outputs.protocol.astype(int)
outputs.visitno= outputs.visitno.astype(float)
outputs.llod= outputs.llod.astype(float)

# MANIFEST COMPLETENESS CHECK
# expecting to see ptid match against visits 4 and 6
manifest3 = pd.read_csv(
    "/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN206/assays/MAAR/misc_files/manifests/hvtn206_batch1_2_reactive_ptid.csv"
)
assert set(outputs.loc[outputs.visitno.isin([4., 6.])].ptid).symmetric_difference(manifest3.ptid) == set()

# expecting ptid-visit combos are an exact match for all visitno==2 in the data
manifest = pd.read_csv(
    "/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN206/assays/MAAR/misc_files/manifests/512-015-2026000119.txt", sep="\t"
)

manifest['pid_vid'] = manifest[['PID','VID']].astype(str).agg("|".join, axis=1)
ptid_visitnos = outputs[['ptid','visitno']].astype(int).astype(str).agg("|".join, axis=1)
assert set(ptid_visitnos[outputs.visitno==2.]).symmetric_difference(manifest.pid_vid) == set()

# SAVE RESULTS SPECIFIC TO THE SEPT UPLOAD
today = datetime.date.today().isoformat()
savedir = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN206/assays/MAAR/misc_files/data_processing/'
outputs.to_csv(savedir + f"MAAR-Results-HVTN206-24Sep2026_SDMC_PROCESSED_{today}.txt", sep='\t', index=False)

# DATA UPLOAD SUMMARIES SPECIFIC TO THE SEPT UPLOAD
# summary of qualitative results
new_summary = outputs.loc[outputs.result_qualitative!='Not Done'].pivot(
    index=['ptid','visitno'],
    columns='assay_subtype',
    values='result_qualitative'
).fillna("Not Done")
new_summary.to_excel(
    savedir + f"HVTN206_MAAR_lab_sept_upload_summary_2026-09-28.xlsx"
)

summary_imputed = outputs.loc[outputs.result_qualitative!='Not Done'].pivot_table(
    index=['ptid','visitno'],
    columns='assay_subtype',
    values='result_qualitative',
    aggfunc='count'
).fillna(0)
summary_imputed.to_excel(
    savedir + f"HVTN206_MAAR_24Sep2026_data_upload_summary_2026-09-29.xlsx"
)

# APPEND PRIOR DATASETS TO THIS NEW ONE
savedir = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN206/assays/MAAR/misc_files/data_processing/'
feb_full_df = pd.read_csv(
    savedir + 'HVTN206_MAAR_full_dataset_processed_2026-02-23.txt',
    sep="\t"
)

assert set(feb_full_df.columns).symmetric_difference(outputs.columns) == set()

new_full = pd.concat([
    feb_full_df,
    outputs
])

new_full.to_csv(savedir + f"HVTN206_MAAR_JanFebSept2026_combined_SDMC_PROCESSED_{today}.txt", sep='\t', index=False)

# CUMULATIVE SAMPLE UPLOAD SUMMARY
cumulative_summary = new_full.loc[new_full.result_qualitative!='Not Done'].pivot(
    index=['ptid','visitno'],
    columns='assay_subtype',
    values='result_qualitative'
).fillna("Not Done")
cumulative_summary.to_excel(
    savedir + f"HVTN206_MAAR_Sept2026_cumulative_sample_summary_{today}.xlsx"
)