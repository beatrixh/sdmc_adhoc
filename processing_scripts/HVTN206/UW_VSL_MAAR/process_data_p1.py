## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 1/19/2026
# Purpose: MAAR
## ---------------------------------------------------------------------------##
import pandas as pd
import numpy as np
import os
import datetime as datetime

import sdmc_tools.process as sdmc_tools
import sdmc_tools.access_ldms as access_ldms

# TODO: clean up notebook formatting

input_data_path="/trials/vaccine/p206/s001/qdata/LabData/MAAR_pass-through/MAAR-Results-07Oct2025-Revised-08Jan2026.xlsx"
df=pd.read_excel(input_data_path)

prev_processed = pd.read_csv(
    "/trials/vaccine/p206/s001/qdata/LabData/MAAR_pass-through/processed_by_sdmc/HVTN206_MAAR_processed_2025-11-14.txt",
    sep="\t"
)

prev_df = pd.read_excel("/trials/vaccine/p206/s001/qdata/LabData/MAAR_pass-through/archive/MAAR-Results-07Oct2025-Final.xlsx")
prev = prev_df['Global Spec IDs'].str.split(",", expand=True).melt(id_vars=None, value_vars=[0,1,2], value_name='guspec')
prev = prev[['guspec']].dropna()
prev.guspec = prev.guspec.str.strip()

set(prev.guspec).difference(df['Global Spec IDs'])

set(df['Global Spec IDs']).difference(prev.guspec)

pd.options.display.max_columns=100

df = df.melt(
    id_vars=['PID','Global Spec IDs'],
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

result_detail = df.result.str.split(",", expand=True).rename(
    columns={0:'result_qualitative', 1:'result_quantitative'}
)
df = pd.concat([
    df.drop(columns='result'),
    result_detail
], axis=1)

# # prev_processed.loc[(
# #     (prev_processed.assay_subtype_lab=="BioRad Geenius HIV-1/2 Ab")
# # )]
# prev_processed.loc[(
#     prev_processed.assay_subtype=="BioRad Geenius HIV-1/2 Ab"
# ), ['result_qualitative','result_quantitative','result_detail','result_units']].drop_duplicates()

#move values in qual column into an 'result_detail' column for biorad geenius
df.loc[(df.assay_subtype_lab=='BioRad Geenius HIV-1/2 Ab') & (df.result_qualitative!='Non-reactive'), 'result_detail'] = df.loc[(df.assay_subtype_lab=='BioRad Geenius HIV-1/2 Ab') & (df.result_qualitative!='Non-reactive'), 'result_qualitative']

# amend biorad geenius qual column to contain "Ab Reactive" for all rows that aren't non-reactive
df.loc[(df.assay_subtype_lab=='BioRad Geenius HIV-1/2 Ab') & (df.result_qualitative.isin(['HIV-1 Ind. gp160','HIV Ind. gp140 and gp160'])), 'result_qualitative'] = 'Ab Reactive'

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
    guspec_col='Global Spec IDs',
    network='hvtn',
    metadata_dict=md,
    ldms=ldms,
)

assert (outputs.ptid.astype(int)!=outputs.pid.astype(int)).sum()==0
outputs.ptid = outputs.ptid.astype(int)

outputs = outputs.drop(
    columns=['global_spec_ids','pid']
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

tmp = outputs.copy()
tmp['guspec_core'] = tmp.guspec.str.rpartition("-")[0]
tmp = tmp.drop(columns='guspec')
tmp = tmp.drop_duplicates()

tmp.shape, prev_processed.shape

tmp = tmp.loc[~(tmp.result_quantitative=="Not Done")]
prev_tmp = prev_processed.loc[~(prev_processed.result_quantitative=="Not Done")]
prev_tmp = prev_tmp.drop(columns='guspec_aliquots_possible')

# tmp.shape, prev_tmp.shape

set(tmp.columns).symmetric_difference(prev_tmp.columns)

sort = tmp.columns.tolist()
prev_tmp = prev_tmp[sort].sort_values(by=sort).reset_index(drop=True)
tmp = tmp[sort].sort_values(by=sort).reset_index(drop=True)

diffs = tmp.compare(prev_tmp)

diffs['instrument_serialno'].drop_duplicates()

today = datetime.date.today().isoformat()
savedir = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN206/assays/MAAR/misc_files/data_processing/'
outputs.to_csv(savedir + f"HVTN206_MAAR_processed_{today}.txt", sep='\t', index=False)

lab_manifest = pd.read_excel(
    '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN206/assays/MAAR/misc_files/lab_submission/lab_submitted_manifest/HVTN206-v2.0.xlsx',
)

# all data guspecs are in lab manifest
set(outputs.guspec).difference(lab_manifest['Global Spec ID']) 
set(lab_manifest['Global Spec ID']).difference(outputs.guspec)


summary = outputs.loc[outputs.result_qualitative!='Not Done']
pivot_summary = pd.pivot_table(
    summary,
    index=['ptid','visitno'],
    columns=['assay_subtype'],
    aggfunc='count',
    values='guspec',
    fill_value=0
)

pivot_summary.to_excel(
    savedir + "HVTN206_sample_summary_2026-01-19.xlsx"
)

summary2 = summary.copy()
summary2['is_reactive'] = False
summary2.loc[summary2.result_qualitative=="Ab Reactive", 'is_reactive'] = True

summary2 = summary2.groupby(['visitno','ptid'])[['is_reactive']].sum().reset_index()
summary2 = summary2.rename(columns={'is_reactive':'count_of_reactives'})
summary2 = summary2.set_index(['visitno','ptid'])

summary2.to_excel(savedir + "HVTN206_MAAR_sample_reactivity_summary_2026-01-19.xlsx")

summary2_og = pd.read_excel(
    '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN206/assays/MAAR/misc_files/data_processing/archive_2026-01-19/HVTN206_MAAR_sample_reactivity_summary.xlsx', index_col=[0,1]
)

summary2_og.reset_index().compare(summary2.reset_index())