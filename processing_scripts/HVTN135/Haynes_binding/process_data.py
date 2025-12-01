## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 11/20/2025
# Purpose: Process the SECOND UPLOAD of ELISA binding data for HVTN 135
## ---------------------------------------------------------------------------##
# CAUTION: the lab initially uploaded HVTN 115 and 135 data from a different run.
# This data was processed by a different script (see HVTN115 dir in this repo).
# This script processes the new upload, appends it to the prior outputs, and saves.
## ---------------------------------------------------------------------------##
import pandas as pd
import numpy as np
import os
import datetime as datetime

import sdmc_tools.process as sdmc_tools
import sdmc_tools.access_ldms as access_ldms

def main():
	# read in new data
	binding_path = "/trials/vaccine/p135/s001/qdata/LabData/elisa_binding_pass-through/HVTN 135 binding master file v3.xlsx"
	binding = pd.read_excel(
	    binding_path,
	    sheet_name='binding summary',
	    skiprows=7
	)

	# these should uniquely identify a row
	assert len(binding) == binding[['GLOBAL_ID','Antigen']].drop_duplicates().shape[0]

	# cast to long on dilution
	binding = binding.melt(
	    id_vars=[
	        'Run Date',
	        'Ag Lot',
	        'Mother-Infant Pair-matched #',
	        'Ag #',
	        'Study ID',
	        'GLOBAL_ID',
	        'Pub ID',
	        'Group',
	        'Treatment',
	        'Subject Type',
	        'PTID',
	        'Visit ID',
	        'Antigen',
	        'Log AUC',
	        'EC50',
	        'End Titer',
	        'Dil50',
	        'Mab Eq'
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

	binding.dilution = binding.dilution.str[1:].astype(int)
	binding['specrole'] = 'Sample'

	binding = binding.rename(columns={
	    'GLOBAL_ID':'guspec',
	    'PTID': 'ptid_lab_submitted'
	})

	# correct typos
	binding.loc[(binding.guspec=='0410-ON1V3E00-001'), 'guspec'] = '0410-0N1V3E00-001'
	binding.loc[(binding.guspec=='0410-ONCQCG00-002'), 'guspec'] = '0410-0NCQCG00-002'
	binding.loc[(binding.guspec=='04100NJV3K00-001'), 'guspec'] = '0410-0NJV3K00-001'

	# pull ldms
	ldms = access_ldms.pull_one_protocol('hvtn', 135)
	ldms = ldms.loc[ldms.guspec.isin(binding.guspec.tolist())].drop_duplicates()

	# standard processing
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

	outputs = sdmc_tools.standard_processing(
	    input_data=binding,
	    input_data_path=binding_path,
	    guspec_col='guspec',
	    network='hvtn',
	    metadata_dict=md,
	    ldms=ldms
	)

	## LDMS CHECKS
	outputs.ptid = outputs.ptid.astype(float)
	outputs.ptid_lab_submitted = outputs.ptid_lab_submitted.astype(float)

	outputs.visitno = outputs.visitno.astype(float)
	outputs.visit_id = outputs.visit_id.astype(float)

	assert (outputs.ptid!=outputs.ptid_lab_submitted).sum()==0
	assert (outputs.visitno!=outputs.visit_id).sum()==0
	assert outputs[['ag_#','antigen']].drop_duplicates().shape[0] == outputs.antigen.nunique()

	outputs = outputs.drop(columns=[
	    'ptid_lab_submitted',
	    'visit_id',
	    'ag_#',
	    'mother-infant_pair-matched_#',
	    'pub_id',
	    'group',
	    'study_id',
	    'treatment',
	    'subject_type',
	    'dil50', #rob confirmed via email to drop these
	    'mab_eq', #rob confirmed via email to drop these
	])

	outputs = outputs.rename(columns={
	    'ag_lot':'analyte_lot',
	    'antigen':'analyte',
	})


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
	    # 'control_name',				#these cols only in first submission
	    # 'control_concentration',
	    # 'control_concentration_units',
	    'dilution',
	    'result',
	    'result_units',
	    'ec50',
	    'end_titer',
	    'log_auc',
	    # 'plate_number',
	    # 'sample_position',
	    'sdmc_processing_datetime',
	    'sdmc_data_receipt_datetime',
	    'input_file_name',
	]
	assert set(outputs.columns).symmetric_difference(reorder) == set()
	outputs = outputs[reorder]


	# read in outputs from Nov 12 upload and combine
	previous = pd.read_csv(
	    '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN135/assays/binding_ELISA/misc_files/data_processing/HVTN135_ELISA_binding_original_upload_processed_2025-11-12.txt',
	    sep="\t"
	)

	outputs.protocol = outputs.protocol.astype(float)
	outputs.run_date = outputs.run_date.astype(str)

	both_outputs = pd.concat([
	    previous,
	    outputs
	])

	new_reorder = [
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
	both_outputs = both_outputs[new_reorder]

	today = datetime.date.today().isoformat()
	savedir135 = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN135/assays/binding_ELISA/misc_files/data_processing/'
	both_outputs.to_csv(savedir135 + f'HVTN135_ELISA_binding_processed_{today}.txt', sep="\t", index=False)


	# read in 115 binding so can create full upload/protocol pivot summary
	outputs115 = pd.read_csv(
	    "/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN115/assays/binding_ELISA/misc_files/data_processing/HVTN115_ELISA_binding_processed_2025-11-12.txt",
	    sep='\t'
	)

	outputs_all = pd.concat([
	    both_outputs, #135
	    outputs115
	])

	outputs_all.plate_number = outputs_all.plate_number.astype(str)

	outputs_all = outputs_all.drop_duplicates() # Standards and Background were saved in both 115 and 135

	# # visual confirm that for one ptid-visit-analyte-upload there are always 12 dilutions
	# outputs_all.loc[outputs_all.specrole=='Sample'].groupby([
	#     'ptid',
	#     'visitno',
	#     'analyte',
	#     'input_file_name',
	# ]).dilution.nunique().value_counts()

	# these are the same
	outputs_all.analyte = outputs_all.analyte.str.strip()
	outputs_all.analyte = outputs_all.analyte.str.replace("Strep","strep")

	samples_summary = pd.pivot_table(
	    outputs_all.loc[(outputs_all.specrole=='Sample') & (outputs_all.dilution==30.)],
	    index=['ptid','visitno'],
	    columns=['analyte'],
	    aggfunc='count',
	    fill_value=0,
	    dropna=False,
	    values='result'
	)

	controls_summary = pd.pivot_table(
	    both_outputs.loc[(both_outputs.specrole!='Sample')],
	    index=['specrole','analyte'],
	    columns=['control_concentration'],
	    aggfunc='count',
	    dropna=False,
	    values='result'
	).dropna(how='all').fillna(0)

	# Define the file path
	savedir115 = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN115/assays/binding_ELISA/misc_files/data_processing/'

	# Create an ExcelWriter object
	with pd.ExcelWriter(savedir135 + "HVTN115_135_ELISA_binding_summary_2025-11-24.xlsx", engine='openpyxl') as writer:
	    # Write each dataframe to a different sheet
	    samples_summary.to_excel(writer, sheet_name='Sample summary', index=True) # index=False prevents writing the DataFrame index to Excel
	    controls_summary.to_excel(writer, sheet_name='Control summary', index=True)

	# # Create an ExcelWriter object
	# with pd.ExcelWriter(savedir115 + "HVTN115_135_ELISA_binding_summary.xlsx", engine='openpyxl') as writer:
	#     # Write each dataframe to a different sheet
	#     samples_summary.to_excel(writer, sheet_name='Sample summary', index=True) # index=False prevents writing the DataFrame index to Excel
	#     controls_summary.to_excel(writer, sheet_name='Control summary', index=True)

if __name__=="__main__":
	main()