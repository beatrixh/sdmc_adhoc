import pandas as pd
import numpy as np
import os
import datetime as datetime

import sdmc_tools.process as sdmc_tools
import sdmc_tools.access_ldms as access_ldms

new_path = "/trials/vaccine/p135/s001/qdata/LabData/elisa_binding_pass-through/251209 HVTN135 Infant Baseline vs p5 Binding.xlsx"
new = pd.read_excel(
    new_path,
    sheet_name='Binding Summary',
    skiprows=12
)

new = new.melt(
    id_vars=[
        'Run Date',
        'Plate ',
        'sample position',
        'Plate Designation',
        'Antigen Lot',
        'Study',
        'Original ID',
        'Subject type',
        'Treatment',
        'Pub ID',
        'Subject ID',
        'Group',
        'Visit',
        'Rep #',
        'Antigen',
        'Log AUC',
    ],
    value_vars=[
        '_30',
        '_90',
        '_270',
        '_810',
        '_2430',
        '_7290',
        '_21870',
        '_65610',
        '_196830',
        '_590490',
        '_1771470',
        '_5314410',
    ],
    value_name='result',
    var_name='dilution',
)

new.dilution = new.dilution.str[1:].astype(int)
new = new.rename(columns={
    'Original ID':'guspec',
    'Subject ID': 'ptid_lab_submitted',
    'Visit':'visitno_submitted',
    'Rep #':'replicate',
    'Antigen':'analyte',
    'Antigen Lot':'analyte_lot',
    'Plate ':'plate_number',
})
new['specrole'] = 'Sample'

md = {
    'network':'HVTN',
    'upload_lab_id': 'BH',
    'assay_lab_name':'Haynes Lab (Duke)',
    'assay_type': 'ELISA',
    'assay_subtype':'Binding',
    'instrument':'SpectraMax',
    'lab_software':'Softmax',
    'lab_software_version':'Softmax 5.3',
    'result_units':'Absorbance (450 nm)',
    'assay_precision':'Quantitative',
}

ldms = access_ldms.pull_one_protocol('hvtn', 135)
ldms = ldms.loc[ldms.guspec.isin(new.guspec.tolist())].drop_duplicates()

outputs_new = sdmc_tools.standard_processing(
    input_data=new,
    input_data_path=new_path,
    guspec_col='guspec',
    network='hvtn',
    metadata_dict=md,
    ldms=ldms
)

# check what this contains -- 
# outputs_new.plate_designation.unique()
#  array(['purple A2', 'green A2', 'purple A1', 'green A1', 'purple A3','green A3'], dtype=object)

# LDMS checks
outputs_new.ptid = outputs_new.ptid.astype(float)
outputs_new.ptid_lab_submitted = outputs_new.ptid_lab_submitted.astype(float)

outputs_new.visitno = outputs_new.visitno.astype(float)
outputs_new.visitno_submitted = outputs_new.visitno_submitted.astype(float)

assert (outputs_new.ptid!=outputs_new.ptid_lab_submitted).sum()==0
assert (outputs_new.visitno!=outputs_new.visitno_submitted).sum()==0


# formatting
outputs_new = outputs_new.drop(columns=[
    'ptid_lab_submitted',
    'visitno_submitted',
    'group',
    'study',
    'treatment',
    'subject_type',
    'plate_designation',
    'pub_id',
])

# # make sure there are 12 dilutions per guspec so we can subset to only one of the dilutions
# outputs_new.groupby(['guspec']).dilution.nunique().value_counts()

# outputs_new.dilution.nunique() == 12

# pivot summary
new_summary = pd.pivot_table(
    data=outputs_new.loc[outputs_new.dilution==30],
    index=['ptid'],
    columns=['analyte','visitno'],
    aggfunc='count',
    fill_value=0,
    values='sdmc_processing_datetime'
)

new_summary.to_excel(
        "/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN135/assays/binding_ELISA/misc_files/data_processing/HVTN135_binding_12_11_upload_summary.xlsx"
)

# read in prior outputs
outputs_old = pd.read_csv(
    '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN135/assays/binding_ELISA/misc_files/data_processing/HVTN135_ELISA_binding_processed_2025-11-25.txt',
    sep='\t'
)

# which columns are we missing in the new one
set(outputs_old.columns).difference(outputs_new.columns)
# {'control_concentration',
#  'control_concentration_units',
#  'control_name',
#  'ec50',
#  'end_titer'}

# which columns are new
set(outputs_new.columns).difference(outputs_old.columns)
# {'replicate'}

outputs = pd.concat([outputs_new, outputs_old])
reorder = [
    'network',
    'protocol',
    'specrole',
    'guspec',
    'ptid',
    'visitno',
    'drawdt',
    'run_date',
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
    'lab_software',
    'lab_software_version',
    'analyte',
    'analyte_lot',
    'control_name',
    'control_concentration',
    'control_concentration_units',
    'dilution',
    'replicate',
    'result',
    'result_units',
    'ec50',
    'end_titer',
    'log_auc',
    'plate_number',
    'sample_position',
    'sdmc_processing_datetime',
    'sdmc_data_receipt_datetime',
    'input_file_name',
]
outputs = outputs[reorder]

outputs.run_date = pd.to_datetime(outputs.run_date).astype(str)
today = datetime.date.today().isoformat()
savedir135 = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN135/assays/binding_ELISA/misc_files/data_processing/'
outputs.to_csv(savedir135 + f'HVTN135_ELISA_binding_processed_{today}.txt', sep="\t", index=False)

infant_list = pd.read_excel("/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN135/assays/HVTN135_infants.xlsx")

# seeing all infant ptids except this one
set(infant_list.PTID.astype(int)).difference(outputs_new.ptid.astype(int))
# {782566697}
