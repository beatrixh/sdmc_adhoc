## ---------------------------------------------------------------------------##
# Author: Sara Thiebaud
# Date: 18 March 2026
# Purpose: HVTN 302 BLI data from Tomaras lab
## ---------------------------------------------------------------------------##
import pandas
import numpy as np
import os
import datetime as datetime
import pdb

import sdmc_tools.process as sdmc_tools
import sdmc_tools.access_ldms as access_ldms

NETWORK = 'hvtn'
PROTOCOL = 302

WORKING_DIR = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN302/assays/BLI/misc_files/data_processing/'
QDATA_FOLDER = '/trials/vaccine/p302/s001/qdata/LabData/BLI_pass-through'
QDATA_FILES = [
    'HVTN_302_BLI_GT_BG505_MD39.3_CD4KO4_pHLsecAvi_20260317.txt',
    'HVTN_302_BLI_GT_BG505_MD39.3_pHLsecAvi_20260317.txt'
]

COLUMN_NAMING = {
    'analyte':'sample_description',
    'analyte_type':'assay_detail',
    'lab':'upload_lab_id',
    'experiment_type':'assay_name',
    'response':'result_quantitative',
    'response_unit':'result_units',
    'response_standard_deviation':'result_stdev',
    'response_coefficient_of_variation':'result_cv',
    'kinetic_off_rate_standard_deviation':'kinetic_off_rate_stdev',
    'kinetic_off_rate_coefficient_of_variation':'kinetic_off_rate_cv',
    'kinetic_dissociation_area_under_curve':'kinetic_dissociation_auc',
    'kinetic_dissociation_area_under_curve_unit':'kinetic_dissociation_auc_units',
    'kinetic_off_rate_unit':'kinetic_off_rate_units',
    'positivity':'result_qualitative',
    'quantifiable':'result_quantifiable',
    'lower_limit_of_quantitation':'lloq',
    'dataset_creation_date':'dataset_creation_date_lab',
    'note':'lab_comments'
}

METADATA = {
    'specrole':'Sample',
    'assay_lab_name':'Tomaras Lab (Duke)',
    'assay_subtype':'Avidity',
    'assay_precision':'Quantitative',
    'instrument':'BLI Octet'
    }

OUTPUT_COLUMNS = [
    'network',
    'protocol',
    'specrole',
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
    'instrument',
    'assay_name',
    'assay_subtype',
    'assay_precision',
    'assay_detail',
    'ligand',
    'ligand_type',
    'sample_description',
    'exclude',
    'experiment_identifier',
    'plate_identifier',
    'result_quantifiable',
    'result_quantitative',
    'result_cv',
    'result_stdev',
    'result_units',
    'kinetic_dissociation_auc',
    'kinetic_dissociation_auc_units',
    'positivity_threshold',
    'result_qualitative',
    'lloq',
    'kinetic_off_rate',
    'kinetic_off_rate_cv',
    'kinetic_off_rate_stdev',
    'kinetic_off_rate_units',
    'lab_comments',
    'value_summary',
    'value_summary_method',
    'value_summary_type',
    'dataset_creation_date_lab',
    'dataset_version',
    'input_file_name',
    'sdmc_processing_datetime',
    'sdmc_data_receipt_datetime',
]

ds = pandas.DataFrame()
ldms = access_ldms.pull_one_protocol(NETWORK, PROTOCOL)

for f in QDATA_FILES:

    fp = os.path.join(QDATA_FOLDER, f)
    df = pandas.read_csv(fp, dtype=str, sep='\t')

    df = df.rename(columns=COLUMN_NAMING)
    df['guspec'] = df.apply(lambda row: row['sample_description'].split(';')[0], axis=1)
    dfp = sdmc_tools.standard_processing(
        input_data=df,
        input_data_path=fp,
        guspec_col='guspec',
        network='hvtn',
        metadata_dict=METADATA,
        ldms=ldms,
    )

    ds = pandas.concat([ds, dfp])

# assert (ds.ptid.astype(int)!=ds.sample_identifier.astype(int)).sum()==0
# assert (ds.visitno.astype(float)!=ds.visit_identifier.astype(float)).sum()==0

outputs = ds.drop(
    columns=['sample_type', 'sample_identifier', 'visit_identifier', 'study']
)

# reorder = [
#     'network',
#     'protocol',
#     'guspec',
#     'specrole',
#     'upload_lab_id',
#     'assay_lab_name',
#     'assay_name_haws',
#     'ptid',
#     'visitno',
#     'drawdt',
#     'spectype',
#     'spec_primary',
#     'spec_additive',
#     'spec_derivative',
#     'assay_name',
#     'assay_subtype',
#     'assay_precision',
#     'lab_software',
#     'instrument',
#     'instrument_serialno',
#     'result_qualitative',
#     'result_quantitative',
#     'result_detail',
#     'result_units',
#     'llod',
#     'sdmc_processing_datetime',
#     'sdmc_data_receipt_datetime',
#     'input_file_name',
# ]

# set(reorder).symmetric_difference(outputs.columns)

# outputs = outputs[reorder]


today = datetime.date.today().isoformat()
outputs.to_csv(WORKING_DIR + f'HVTN302_BLI_processed_{today}.txt', sep='\t', index=False)

# outputs[['result_qualitative','result_detail','result_units']].drop_duplicates()

# outputs[[i for i in outputs.columns if 'result' in i]].drop_duplicates()
summary = pandas.pivot_table(
    outputs,
    index=['ptid','visitno'],
    columns='ligand',
    aggfunc='count',
    values='input_file_name',
    dropna=False
).dropna(how='all').fillna(0)

# summary.to_excel(os.path.join(WORKING_DIR, f'HVTN302_BLI_summary_{today}.xlsx'))


manifest = pandas.read_excel(os.path.join(WORKING_DIR, 'manifests', 'HVTN 302 v8 and v12 Serum Shipped to the Tomaras Lab as of 23Mar2026.xlsx'))

manifest = manifest.rename(columns={'Original ID':'guspec'})
m = pandas.merge(outputs, manifest, how='outer', on='guspec', indicator=True)
# m.to_csv('manifest_comparison_' + today + '.csv', index=False)
