## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 03/02/2026
# Purpose: HVTN142 RNA Sequencing - IQVIA / Vir
## ---------------------------------------------------------------------------##
import pandas as pd
import numpy as np
import os
import datetime as datetime

import sdmc_tools.process as sdmc_tools
import sdmc_tools.access_ldms as access_ldms

fpath = '/trials/vaccine/p142/s001/qdata/LabData/RNASeq_pass-through/2024-06-27_1388-V101_34663/Samples.2024-06-27_1388-V101_34663_23912363_RNA01.contents.txt'

data = pd.read_csv(fpath)

tmp = []
for i, row in data.iterrows():
    tmp += [row.values[0].strip().split(" ")]

data = pd.DataFrame(tmp, columns=['file_size_bytes','assay_date','hex_string','filepath'])

folder_dir = '/trials/vaccine/p142/s001/qdata/LabData/RNASeq_pass-through/2024-06-27_1388-V101_34663/Samples.2024-06-27_1388-V101_34663_23912363_RNA01/'
files = []
for d in os.listdir(folder_dir):
    for f in os.listdir(folder_dir + d):
        files += [d + "/" + f]

print(set(data.filepath).difference(files))

# there's this one row that doesn't seem to be in the manifest
print(set(files).difference(data.filepath))

# check that file_size_bytes is what we expected
sizes = [os.path.getsize(folder_dir + i) for i in data.filepath]
(pd.Series(sizes).astype(float) - data.file_size_bytes.astype(float)).max() == 0.0

# merge on row of missing data
new_row = pd.DataFrame({
    'file_size_bytes':[os.path.getsize(folder_dir + '20240806-24-2110707.FASTQs/0373-01T7ZC00-001_L1.UDI_0001.adapters')],
    'filepath':['20240806-24-2110707.FASTQs/0373-01T7ZC00-001_L1.UDI_0001.adapters'],
    'sdmc_qc_flag':[True],
    'sdmc_qc_flag_detail':["This file wasn't included in this manifest by the lab. We've merged it on because the lab did submit this file."]
})

data['sdmc_qc_flag'] = False

data = pd.concat([
    data,
    new_row
])

data = pd.concat([
    data,
    data.filepath.str.split("/", expand=True).rename(columns={0:'folder_name',1:'filename'})
], axis=1)

data['guspec'] = data.filename.str.split(".", expand=True)[0].str[:-3]
data['file_extension'] = data.filename.str[14:].str.partition('.')[2].str.partition('.')[2]
data['udi_id'] = data.filename.str[21:29]
data['read_number'] = data.filename.str[27:].str.partition(".")[0].str.partition("_")[2]

# # check how these are distributed
# data.groupby('file_extension').count()

# set(data.loc[(data.guspec=='0373-01VMQD00-001')].row_identifier).difference(data.loc[(data.guspec=='0373-01T7ZC00-001')].row_identifier)


## pull ldms, standard processing
ldms = access_ldms.pull_one_protocol('hvtn', 142)

md = {
    'network': 'HVTN',
    'protocol': 142,
    'upload_lab_id': 'N/A',
    'assay_lab_name': 'Q2 (IQVIA)',
    'assay_type': 'RNA Sequencing',
    'assay_subtype': 'mRNA',
    'assay_kit': 'Illumina TruSeq',
    'instrument': 'Illumina NovaSeq 6000',
    # 'lab_software': 'SOFTWARE',
    # 'lab_software_version': 'SOFTWARE VERSION',
    'assay_precision': 'Quantitative',
}

outputs = sdmc_tools.standard_processing(
    input_data=data,
    input_data_path=fpath,
    guspec_col="guspec",
    network="HVTN",
    metadata_dict=md,
    ldms=ldms
)

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
    'assay_kit',
    'instrument',
    'assay_date',
    'assay_precision',
    'filepath',
    'folder_name',
    'filename',
    'file_extension',
    'read_number',
    'udi_id',
    'hex_string',
    'file_size_bytes',
    'sdmc_processing_datetime',
    'sdmc_data_receipt_datetime',
    'input_file_name',
    'sdmc_qc_flag',
    'sdmc_qc_flag_detail',
]
assert set(outputs.columns).symmetric_difference(reorder) == set()
outputs = outputs[reorder]

savedir = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN142/assays/RNASeq/misc_files/data_processing/'
today = datetime.date.today().isoformat()
outputs.to_csv(savedir + f"DRAFT_HVTN142_RNASeq_processed_{today}.txt", sep="\t", index=False)

# pivot summary
sample_summary = pd.pivot_table(
    outputs,
    index=['ptid'],
    columns=['visitno'],
    values='file_size_bytes',
    aggfunc='count',
    fill_value=0
)

today = datetime.date.today().isoformat()
sample_summary.to_excel(
    savedir + f"HVTN142_RNASeq_file_summary_{today}.xlsx"
)

## manifest check
manifest1 = pd.read_csv(
    '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN142/assays/RNASeq/misc_files/manifests/481-999097-0000057255.txt',
    sep='\t'
)

manifest2 = pd.read_csv(
    '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN142/assays/RNASeq/misc_files/manifests/512--2024000130.txt',
    sep='\t'
)

set(manifest1.GLOBAL_ID).difference(outputs.guspec)

set(manifest2.GLOBAL_ID).difference(outputs.guspec)

set(outputs.guspec).difference(set(manifest1.GLOBAL_ID).union(manifest2.GLOBAL_ID))

manifest2.loc[manifest2.GLOBAL_ID.isin(['0374-01GLMG00-001', '0378-01064A00-001', '0378-0106LJ00-001'])]

manifest2.loc[manifest2.GLOBAL_ID.isin(['0374-01GLMG00-001', '0378-01064A00-001', '0378-0106LJ00-001']),
['GLOBAL_ID','PID','VID']
]