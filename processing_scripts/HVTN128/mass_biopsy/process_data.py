import pandas as pd
import numpy as np
import os
import sdmc_tools.access_ldms as access_ldms
import sdmc_tools.process as sdmc_tools
import datetime

input_data_path = "/trials/vaccine/p128/s001/qdata/LabData/Specimen_Mass_pass-through/HVTN_128_Biopsy_Mass_GT_20251120_B.xlsx"
data = pd.read_excel(
    input_data_path,
)

# # make sure lab, network, study, mass_units all constant
# data.nunique()


# all data guspecs are in ldms (128)
ldms = access_ldms.pull_one_protocol('hvtn', 128)
set(data.global_identifier).difference(ldms.guspec)
ldms = ldms.loc[ldms.guspec.isin(data.global_identifier)].drop_duplicates()

# LDMS CHECKS
compare = data.merge(
    ldms,
    left_on='global_identifier',
    right_on='guspec',
    how='left',
)

assert len(compare) == len(data)

assert (compare.txtpid.astype(int)!=compare.sample_identifier.astype(int)).sum() == 0
assert (compare.visit_identifier.astype(int)!=compare.vidval.astype(int)).sum() == 0
assert (compare.lstudy.astype(int)!=compare.study.astype(int)).sum() == 0

compare['sample_type_ldms'] = compare.primstr+","+compare.addstr+","+compare.dervstr
assert (compare.sample_type_ldms!=compare.sample_type).sum() == 0

# drop these, merge them on from LDMS
data = data.drop(columns=[
    'lab','network','sample_identifier','sample_type','study','visit_identifier'
])
data = data.rename(columns={
    'global_identifier':'guspec',
    'mass':'result',
    'mass_units':'result_units',
})

# standard formatting
md = {
    'network':'HVTN',
    'upload_lab_id':'GT',
    'assay_lab_name':'Tomaras Lab (Duke)',
    'assay_type':'Specimen Metadata',
    'assay_subtype':'Biopsy Mass',
    'assay_precision':'Quantitative',
    'specrole':'Sample',
}

outputs = sdmc_tools.standard_processing(
    input_data=data,
    input_data_path=input_data_path,
    guspec_col='guspec',
    network='hvtn',
    metadata_dict=md,
    ldms=ldms
)

len(data) == len(outputs)

# save to .txt
today = datetime.date.today().isoformat()
savedir = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN128/assays/mass_biopsy/misc_files/data_processing/'
outputs.to_csv(savedir + f'HVTN128_mass_biopsy_processed_{today}.txt', index=False, sep="\t")