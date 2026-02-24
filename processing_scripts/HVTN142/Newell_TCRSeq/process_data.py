## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 23 Feb 2026
# Purpose: Process HVTN142 TCRSeq metadata file for qdata
## ---------------------------------------------------------------------------##

# TODO -- clean up and update once we've finalized format

import pandas as pd
import numpy as np
import os
import datetime
import sdmc_tools.process as sdmc
import sdmc_tools.access_ldms as access_ldms

# input_data_path = "/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN142/assays/TCRSequencing/misc_files/HVTN_142_sample_prep_infomation.xlsx"
input_data_path = '/trials/vaccine/p142/s001/qdata/LabData/TCR_pass-through/20260223-01/HVTN_142_sample_prep_infomation.xlsx'

df = pd.read_excel(input_data_path)

datadir = '/trials/vaccine/p142/s001/qdata/LabData/TCR_pass-through/'
fnames = os.listdir(datadir)

fnames = list(set(fnames).difference(['archive','20260116-03','20260223-01']))
fnames = pd.DataFrame({'fname':fnames})
fnames = pd.concat([fnames, fnames.fname.str.split(".", expand=True)], axis=1).iloc[:,:3]
fnames.columns = ['mixcr_file_name', 'guspec', 'locus']

fnames.locus = fnames.locus.map({
    'clones_TRG': 'TRG',
    'clones_TRA': 'TRA',
    'clones_TRB': 'TRB',
    'clones_TRD': 'TRD',
    'clones_': 'N/A'
})
set(fnames.guspec).difference(df.guspec)

set(fnames.guspec).difference(df.mixcr_file_name)

test = df.drop(columns='mixcr_file_name').merge(
    fnames,on='guspec',how='left'
)

test.shape, df.shape

fnames.groupby(['guspec']).nunique().value_counts()

fnames.guspec.nunique(), df.guspec.nunique()

set(fnames.guspec).difference(df.guspec)

82 * 4 + 1 + (144-82)

df.mixcr_file_name.isna().sum(), df.mixcr_file_name.notna().sum()

# there's one missing results, so 62 and 82
62 + (82)*4 + 1

df = df.drop(columns='mixcr_file_name').merge(
    fnames,on='guspec',how='left'
)

# only seeing Dec data. are we expecting Oct data?
df.Sequencing_date.value_counts(dropna=False)

run_id_map = {
    '2025-12-04':'251204_VH00699_598_AACWY3THV-1',
    'Oct': '251028_VH00738_347_AAHNFKMM5',
}

df['lab_run_id'] = df.Sequencing_date.astype(str).map(run_id_map)

df = df.rename(columns={
    'templates':'total_reads_per_guspec',
})

df.loc[df.has_result.isna(),'has_result'] = 'NA'
df['run_outcome'] = df.has_result.map({
    'Y':'Sequencing succeeded',
    'N':'Sequencing failed',
    'NA':'Not sequenced',
})

ldms = access_ldms.pull_one_protocol('hvtn', 142)

md = {
    'network': 'HVTN',
    'protocol': 142,
    'upload_lab_id': 'N9',
    'assay_lab_name': 'Newell Lab (Fred Hutch)',
    'assay_type': 'TCR Sequencing',
    'assay_subtype': 'NGS',
    'assay_kit_name':'Cellecta DriverMap AIR TCR-BCR Profiling Kit (Human DNA)',
    'assay_kit_catalog_number':'DMAIR2-HTD-96',
    'instrument': 'Illumina NextSeq 2000',
    'lab_software': 'MiXCR',
    # 'lab_software_version': 'MiXCR 4.0 or later',
    'assay_precision': 'Quantitative',
    'total_dna_input_mass': '2 micrograms'
}

outputs = sdmc.standard_processing(
    input_data=df,
    input_data_path=input_data_path,
    guspec_col="guspec",
    network="HVTN",
    metadata_dict=md,
    ldms=ldms
)

# we expect these to be missing, michelle said they were excluded
excluded_subjects = [
    761629758,
    778239504,
    787322302,
    821107655,
    821447932,
    863535806,
]

outputs.loc[outputs.ptid.astype(int).isin(excluded_subjects),'run_outcome_detail'] = "Subject ID had a greater than 50pct failure during library prep. Not taken to sequencing stage"
outputs.loc[~(outputs.ptid.astype(int).isin(excluded_subjects)) & (outputs.sequence_file_name.isna()),'run_outcome_detail'] = "1-3 libraries did not meet sequencing standards, excluded fin the library pool for sequencing"
outputs.loc[outputs.has_result=="N", 'run_outcome_detail'] = outputs.loc[outputs.has_result=="N", 'comments']

outputs = outputs.rename(columns={
    'lab_data_processing_version':'lab_software_version'
})

Convergence failed in the estimation of the normalization factor for this sample during post-mixcr template estimation, and thus could not be collapsed against reads. This sample should be considered failed.

failed_guspec = '0169-02DPKK00-001'
failed_missingness_reason = "Convergence failed in the estimation of the normalization factor for this sample during post-mixcr template estimation, and thus could not be collapsed against reads. This sample should be considered failed."

outputs.loc[outputs.guspec==failed_guspec,'run_outcome_detail'] = failed_missingness_reason
outputs.loc[outputs.guspec==failed_guspec,'comments'] = failed_missingness_reason
outputs.loc[outputs.guspec==failed_guspec,'run_outcome'] = 'Sequencing failed'

outputs.loc[outputs.mixcr_file_name=="0378-010C3E00-001.clones_nan.quantified.tsv", 'run_outcome'] = "N/A"
outputs.loc[outputs.mixcr_file_name=="0378-010C3E00-001.clones_nan.quantified.tsv", 'run_outcome_detail'] = "This is an empty file; the lab said we could disregard it."

outputs.loc[outputs.has_result=="NA",['has_result','run_outcome','run_outcome_detail']].drop_duplicates()

outputs[['has_result','run_outcome','run_outcome_detail']].drop_duplicates()

outputs = outputs.drop(columns=[
    'sequence_file_name',
    'sequencing_method',
    'comments',
])

reorder = [
    'network',
    'protocol',
    'guspec',
    'ptid',
    'visitno',
    'drawdt',
    'spectype',
    'spec_primary',
    'spec_additive',
    'spec_derivative',
    'upload_lab_id',
    'assay_lab_name',
    'assay_type',
    'assay_subtype',
    'assay_precision',
    'assay_kit_name',
    'assay_kit_catalog_number',
    'instrument',
    'lab_software',
    'lab_software_version',
    'total_dna_input_mass',
    'locus',
    'sequencing_date',
    'has_result',
    'run_outcome',
    'run_outcome_detail',
    'index_id',
    'index_sequence_1',
    'index_sequence_2',
    'mixcr_file_name',
    'sample_preparation',
    'sequence_type',
    'genomic_region',
    # 'sequence_file_name',
    # 'sequencing_method',
    'total_reads_per_guspec',
    # 'comments',
    'lab_run_id',
    'sdmc_processing_datetime',
    'sdmc_data_receipt_datetime',
    'input_file_name',
]
set(outputs.columns).symmetric_difference(reorder)

assert set(outputs.columns).symmetric_difference(reorder) == set()
outputs = outputs[reorder]

today = datetime.date.today().isoformat()
savedir = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN142/assays/TCRSequencing/misc_files/data_processing/'
outputs.to_csv(savedir + f"HVTN142_TCRSequencing_metadata_file_processed_{today}.txt", sep="\t", index=False)

prev = pd.read_csv(
    '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN142/assays/TCRSequencing/misc_files/data_processing/HVTN142_TCRSequencing_metadata_file_processed_2026-02-19.txt',
    sep='\t'
)

prev.shape, outputs.shape

for col in ['ptid','visitno']:
    prev[col] = prev[col].astype(str)
    outputs[col] = outputs[col].astype(str)

compare = prev.compare(outputs.drop(columns=['run_outcome']).rename(columns={'run_outcome_detail':'missingness_reason'}))

compare['has_result'].dropna()

prev.iloc[compare['has_result'].dropna().index][['has_result','missingness_reason']]

outputs[['has_result','missingness_reason']].drop_duplicates()

manifest = pd.read_csv(
    '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN142/assays/TCRSequencing/misc_files/manifests/512-373-2025000069.txt', sep="\t"
)

set(manifest.GLOBAL_ID).difference(outputs.guspec)

len(set(outputs.guspec).difference(manifest.GLOBAL_ID))

today = datetime.date.today().isoformat()
summary = pd.pivot_table(
    outputs,
    index='ptid',
    columns='visitno',
    values='mixcr_file_name',
    aggfunc='count',
    fill_value=0
)
summary.to_excel(savedir + f"HVTN142_TCRSeq_metadata_submission_counts_{today}.xlsx")

addmandir = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN142/assays/TCRSequencing/misc_files/manifests/additional_manifests/'
barc_manifest_fname = 'Copy of BARC HVTN 142 PBMC_DNA PAXgene_Tempus for Newell Lab 06May2025.xlsx'
precision_manifest_fname = 'Precision HVTN 142 PBMC_DNA PAXgene for Newell Lab via FHCC 16Apr2025.xlsx'

barc = pd.read_excel(addmandir + barc_manifest_fname)
prec = pd.read_excel(addmandir + precision_manifest_fname)

set(barc['Global Unique Id']).difference(outputs.guspec)

set(prec['Original ID']).difference(outputs.guspec)

set(outputs.guspec).difference(set(barc['Global Unique Id']).union(prec['Original ID']))

## check against projected visits

projvisits = pd.read_csv(
    '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN142/assays/HVTN142_projected_visits_all.tsv',
    sep='\t'
)
projvisits = projvisits.loc[projvisits.Part=="A"]
# set(outputs.ptid.astype(int)).difference(projvisits.PTID.astype(int))
# len(set(projvisits.PTID.astype(int)).difference(outputs.ptid.astype(int)))
visits_all = outputs.visitno.unique().tolist()
# visits_all

# pd.pivot_table(
#     outputs,
#     index=['ptid'],
#     columns='visitno',
#     values='mixcr_file_name',
#     aggfunc='count'
# )

projvisits = projvisits.melt(
    id_vars=['PTID','Enrollment Date','Part','Group']
)

projvisits = pd.concat([projvisits, projvisits.variable.str.split(" ", expand=True).rename(columns={0:'coltype',1:'visitno'})], axis=1)
projvisits= pd.pivot(projvisits, index=['PTID','Enrollment Date', 'Part', 'Group','visitno'],
        columns='coltype', values='value').reset_index()

expected = projvisits.loc[projvisits.Note=='Complete']

expected.visitno = expected.visitno.astype(int)
outputs.visitno = outputs.visitno.astype(int)
expected.PTID = expected.PTID.astype(int)
outputs.ptid = outputs.ptid.astype(int)

expected = expected.loc[expected.visitno.isin([1,  7, 12, 15, 16, 18])]
compare = outputs[['ptid','visitno','drawdt']].merge(expected, left_on=['ptid','visitno'], right_on=['PTID','visitno'], how='outer')
compare.loc[compare.drawdt.isna()]

