## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 6 January 2026
# Purpose: Process HVTN142 TCRSeq metadata file for qdata
## ---------------------------------------------------------------------------##

# TODO -- clean up and update once we've finalized format


import pandas as pd
import numpy as np
import os
import datetime
import sdmc_tools.process as sdmc
import sdmc_tools.access_ldms as access_ldms

input_data_path = "/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN142/assays/TCRSequencing/misc_files/HVTN_142_sample_prep_infomation.xlsx"

df = pd.read_excel(input_data_path)

datadir = '/trials/vaccine/p142/s001/qdata/LabData/TCR_pass-through/'
fnames = os.listdir(datadir)

fnames = pd.DataFrame({'fname':fnames})
fnames = pd.concat([fnames, fnames.fname.str.split(".", expand=True)], axis=1)
fnames.columns = ['mixcr_file_name', 'guspec', 'locus', 'ext']

fnames.locus = fnames.locus.map({
    'clones_TRG': 'G',
    'clones_TRA': 'A',
    'clones_TRB': 'B',
    'clones_TRD': 'D',
    'clones_': 'N/A'
})

fnames = fnames.iloc[:,:-1]

# make sure this merge does what i want
# set(fnames.guspec).difference(df.guspec)

# set(fnames.guspec).difference(df.mixcr_file_name)

# test = df.drop(columns='mixcr_file_name').merge(
#     fnames,on='guspec',how='left'
# )

# test.shape, df.shape

# df.mixcr_file_name.isna().sum(), df.mixcr_file_name.notna().sum()

# # there's one missing results, so 62 and 82
# 62 + (82)*4 + 1

df = df.drop(columns='mixcr_file_name').merge(
    fnames,on='guspec',how='left'
)

ldms = access_ldms.pull_one_protocol('hvtn', 142)

md = {
    'network': 'HVTN',
    'protocol': 142,
    'upload_lab_id': 'N9',
    'assay_lab_name': 'Newell Lab (Fred Hutch)',
    'assay_type': 'TCR Sequencing',
    'assay_subtype': 'NGS',
    'instrument': 'Illumina NextSeq 2000',
    'lab_software': 'MiXCR',
    # 'lab_software_version': 'MiXCR 4.0 or later',
    'assay_precision': 'Qualitative',
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

outputs.loc[outputs.ptid.astype(int).isin(excluded_subjects),'missingness_reason'] = "Subject ID had a greater than 50pct failure during library prep. Not taken to sequencing stage"
outputs.loc[~(outputs.ptid.astype(int).isin(excluded_subjects)) & (outputs.sequence_file_name.isna()),'missingness_reason'] = "1-3 libraries did not meet sequencing standards, excluded fin the library pool for sequencing"
outputs.loc[outputs.has_result=="N", 'missingness_reason'] = outputs.loc[outputs.has_result=="N", 'comments']

outputs = outputs.rename(columns={
    'lab_data_processing_version':'lab_software_version'
})

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
    'instrument',
    'lab_software_version',
    'lab_software',
    'locus',
    'sequencing_date',
    'has_result',
    'missingness_reason',
    'index_id',
    'index_sequence_1',
    'index_sequence_2',
    'mixcr_file_name',
    'sample_preparation',
    'sequence_type',
    'genomic_region',
    # 'sequence_file_name',
    # 'sequencing_method',
    'templates',
    # 'comments',
    'sdmc_processing_datetime',
    'sdmc_data_receipt_datetime',
    'input_file_name',
]

assert set(outputs.columns).symmetric_difference(reorder) == set()
outputs = outputs[reorder]

today = datetime.date.today().isoformat()
savedir = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN142/assays/TCRSequencing/misc_files/data_processing/'
outputs.to_csv(savedir + f"HVTN142_TCRSequencing_metadata_file_processed_{today}.txt", sep="\t", index=False)