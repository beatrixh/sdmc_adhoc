## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 11/24/2025
# Purpose: Process the SECOND UPLOAD of ELISA blocking data for HVTN 135
## ---------------------------------------------------------------------------##
# CAUTION: the lab initially uploaded HVTN 115 and 135 data from a different run.
# This data was processed by a different script (see HVTN115 dir in this repo).
# This script processes the new upload, appends it to the prior outputs, and saves.
## ---------------------------------------------------------------------------##

## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 11/10/2025
# Purpose: Blocking
## ---------------------------------------------------------------------------##
import pandas as pd
import numpy as np
import os
import datetime as datetime

import sdmc_tools.process as sdmc_tools
import sdmc_tools.access_ldms as access_ldms

def main():
	blocking_path = "/trials/vaccine/p135/s001/qdata/LabData/elisa_blocking_pass-through/HVTN 135 blocking master file v3.xlsx"
	blocking = pd.read_excel(
	    blocking_path,
	    sheet_name="summary v2",
	    skiprows=9
	)

	blocking_header = pd.read_excel(
	    blocking_path,
	    sheet_name="summary v2",
	    skiprows=7,
	    nrows=1
	)

	blocking_header = blocking_header.iloc[:,9:]

	targets = [i.strip() for i in blocking_header.columns if 'Unnamed' not in i]
	analytes = blocking_header.dropna(axis=1).iloc[0].tolist()

	b0 = blocking.columns.tolist()
	b00, b01 = b0[:9], b0[9:]
	b01 = [a+"|"+m for (a,m) in zip([i for i in analytes for _ in range(2)], b01)]
	blocking.columns = b00 + b01

	blocking = blocking.melt(
	    id_vars=b00,
	    value_vars=b01,
	    value_name='result',
	    var_name='analyte|metric'
	)

	target_map = {
	    'CH103_UCA': targets[0],
	    'sCD4': targets[0],
	    'CH106': targets[0],
	    'CH235.12': targets[1],
	    'CH235_UCA': targets[1],
	    '2G12': targets[2],
	    'CH01': targets[2],
	    'PGT145': targets[2],
	    'DH270.6': targets[3],
	}

	# split out antigen and metric
	antigen_metrics = blocking['analyte|metric'].str.split("|", expand=True).rename(columns={0:'analyte',1:'metric'})
	blocking = pd.concat([blocking, antigen_metrics], axis=1).drop(columns='analyte|metric')
	blocking.metric = blocking.metric.str.split(".", expand=True)[0]

	blocking["target"] = blocking.analyte.map(target_map)

	blocking= blocking.pivot(
	    columns='metric',
	    index=[
	        'Study ID',
	        'Mother-Infant Pair-matched #',
	        'Global ID',
	        'Group',
	        'Treatment',
	        'Pub ID',
	        'Subject Type',
	        'PTID',
	        'Visit',
	        'target',
	        'analyte',
	    ],
	    values='result'
	).reset_index()

	blocking[['PTID','Pub ID']].drop_duplicates().shape

	blocking = blocking.rename(columns={
	    'Global ID': 'guspec',
	    'PTID': 'ptid_lab_submitted',
	    '% Blocking': 'result_pct_blocking',
	    'Stdev': 'result_pct_blocking_stddev',
	})

	blocking.loc[(blocking.guspec=='0410-ON1V3E00-001'), 'guspec'] = '0410-0N1V3E00-001'
	blocking.loc[(blocking.guspec=='0410-ONCQCG00-002'), 'guspec'] = '0410-0NCQCG00-002'
	blocking.loc[(blocking.guspec=='04100NJV3K00-001'), 'guspec'] = '0410-0NJV3K00-001'

	# drop blanks
	blocking = blocking.loc[(blocking.result_pct_blocking.notna()) | (blocking.result_pct_blocking_stddev.notna())]

	ldms = access_ldms.pull_one_protocol('hvtn', 135)
	assert set(blocking.guspec).difference(ldms.guspec) == set()
	ldms = ldms.loc[ldms.guspec.isin(blocking.guspec.tolist())]

	# standard processing
	md = {
	    'network':'HVTN',
	    'upload_lab_id': 'BH',
	    'assay_lab_name':'Haynes Lab (Duke)',
	    'specrole':'Sample',
	    'assay_type': 'ELISA',
	    'assay_subtype':'Blocking',
	    'instrument':'SpectraMax',
	    'lab_software':'Softmax',
	    'lab_software_version':'Softmax 5.3',
	    'assay_precision':'Quantitative',
	    'result_units':'Absorbance (450 nm)',
	    'control_concentration_units': 'ug/ml'
	}

	outputs = sdmc_tools.standard_processing(
	    input_data=blocking,
	    input_data_path=blocking_path,
	    guspec_col='guspec',
	    network='hvtn',
	    metadata_dict=md,
	    ldms=ldms
	)

	assert (outputs.visit.astype(float)!=outputs.visitno.astype(float)).sum() == 0
	assert (outputs.ptid_lab_submitted.astype(float)!=outputs.ptid.astype(float)).sum() == 0
	assert outputs[['study_id','protocol']].drop_duplicates().shape[0] == 1

	outputs = outputs.drop(columns=['study_id','visit','ptid_lab_submitted'])

	#lab shouldnt be the source for this info
	outputs = outputs.drop(columns=[
	    'group',
	    'mother-infant_pair-matched_#',
	    'pub_id',
	    'subject_type',
	    'treatment',
	])

	reorder = [
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
	    'assay_type',
	    'assay_subtype',
	    'assay_precision',
	    'instrument',
	    'lab_software',
	    'lab_software_version',
	    'analyte',
	    'target',
	    'control_concentration_units',
	    'result_units',
	    'result_pct_blocking',
	    'result_pct_blocking_stddev',
	    'sdmc_processing_datetime',
	    'sdmc_data_receipt_datetime',
	    'input_file_name'
	]

	assert set(reorder).symmetric_difference(outputs.columns) == set()
	outputs = outputs[reorder]

	# processed version of data from lab's initial 135 blocking upload
	previous = pd.read_csv(
	    "/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN135/assays/blocking_ELISA/misc_files/data_processing/HVTN135_ELISA_bnab_blocking_processed_2025-11-12.txt",
	    sep="\t"
	)

	outputs.protocol = outputs.protocol.astype(float)
	outputs.ptid = outputs.ptid.astype(float)
	outputs.visitno = outputs.visitno.astype(float)

	final = pd.concat([previous, outputs])

	savedir = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN135/assays/blocking_ELISA/misc_files/data_processing/'
	today = datetime.date.today().isoformat()

	# save new 135 blocking outputs
	final.to_csv(
	    savedir + f"HVTN135_ELISA_blocking_processed_{today}.txt",
	    sep='\t',
	    index=False
	)


	# generate summary across uploads and protocols -------------------- ##

	outputs115 = pd.read_csv(
	    "/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN115/assays/blocking_ELISA/misc_files/data_processing/HVTN115_ELISA_bnab_blocking_processed_2025-11-12.txt",
	    sep='\t'
	)

	outputs_all = pd.concat([
	    final,
	    outputs115
	])

	assert len(outputs_all.loc[outputs_all.specrole=='Sample']) == len(outputs_all.loc[outputs_all.specrole=='Sample'].drop_duplicates())
	assert len(outputs_all.loc[outputs_all.specrole!='Sample']) , 2*len(outputs_all.loc[outputs_all.specrole!='Sample'].drop_duplicates())

	# get rid of duplicates standards (saved both to 115 and 135)
	outputs_all = outputs_all.drop_duplicates()

	# make sure guspec/input file/analyte uniquely identify rows
	outputs_all.loc[outputs_all.specrole=='Sample'].groupby(['guspec','input_file_name','analyte']).count().sdmc_data_receipt_datetime.value_counts()

	# visual check these arent dups of one another
	outputs_all.analyte.unique().tolist()

	samples_summary = pd.pivot_table(
	    outputs_all.loc[(outputs_all.specrole=='Sample')],
	    index=['ptid','visitno'],
	    columns=['analyte'],
	    aggfunc='count',
	    fill_value=0,
	    dropna=False,
	    values='result_pct_blocking'
	)

	# make sure each row accounted for
	assert samples_summary.sum().sum() == len(outputs_all.loc[(outputs_all.specrole=='Sample')])

	controls_summary = pd.pivot_table(
	    outputs_all.loc[(outputs_all.specrole!='Sample')],
	    index=['specrole','analyte'],
	    columns=['control_concentration'],
	    aggfunc='count',
	    dropna=False,
	    values='result'
	).dropna(how='all').fillna(0)


	# Create an ExcelWriter object
	with pd.ExcelWriter(savedir + "HVTN115_135_ELISA_blocking_summary_2025-11-24.xlsx", engine='openpyxl') as writer:
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