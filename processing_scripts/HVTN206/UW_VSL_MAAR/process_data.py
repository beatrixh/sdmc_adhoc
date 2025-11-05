## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 10/22/2025
# Purpose: MAAR ad hoc processing
# Note: originally code split guspec into columns 'guspec1', 'guspec2', 
# 'guspec3'; updated code to split into 'guspec_core' and 'guspec_aliquots' instead.
# code still starts off with the 1/2/3 columns.
## ---------------------------------------------------------------------------##
import pandas as pd
import numpy as np
import os
import datetime as datetime

import sdmc_tools.process as sdmc_tools
import sdmc_tools.access_ldms as access_ldms

def main():
	# read in data
	input_path = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN206/assays/MAAR/misc_files/lab_submission/MAAR-Results-07Oct2025-Final.xlsx'
	df = pd.read_excel(input_path)

	# cast to long along assay subtype
	df = df.melt(
	    id_vars=['PID','Global Spec IDs'],
	    value_vars=['Abbott HIV-1/2 Ag/Ab',
	       'BioRad HIV-1/2 Ag/Ab Combo', 'BioRad HIV-1/2 Ab +O (3rd generation)',
	       'Alere Determine HIV-1/2 Ag/Ab', 'Diasorin Liaison HIV-1/2 Ag/Ab',
	       'BioRad Geenius HIV-1/2 Ab', 'INSTI HIV-1/2 Ab Rapid',
	       'OraQuick HIV-1/2 Ab Rapid', 'Abbott HIV-1 RNA'],
	    var_name='assay_subtype_lab',
	    value_name='result'
	)

	# split concatenated guspec lists into their own columns
	guspec_cols = df['Global Spec IDs'].str.replace(" ","").str.split(",", expand=True)
	guspec_cols.columns = ['guspec1','guspec2','guspec3']
	df = pd.concat([guspec_cols, df.drop(columns='Global Spec IDs')], axis=1)

	# split result column out into qual/quant
	df = pd.concat([df.drop(columns='result'), df.result.str.split(",", expand=True).rename(columns={0:'result_qualitative', 1:'result_quantitative'})], axis=1)

	#move values in qual column into an 'result_detail' column for biorad geenius
	df.loc[(df.assay_subtype_lab=='BioRad Geenius HIV-1/2 Ab') & (df.result_qualitative!='Non-reactive'), 'result_detail'] = df.loc[(df.assay_subtype_lab=='BioRad Geenius HIV-1/2 Ab') & (df.result_qualitative!='Non-reactive'), 'result_qualitative']

	# amend biorad geenius qual column to contain "Ab Reactive" for all rows that aren't non-reactive
	df.loc[(df.assay_subtype_lab=='BioRad Geenius HIV-1/2 Ab') & (df.result_qualitative!='Non-reactive'), 'result_qualitative'] = 'Ab Reactive'

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

	# df.loc[(df.assay_subtype_lab=='Abbott HIV-1 RNA'), 'llod'] = '40'

	# merge on metadata from Sara
	maar_metadata = pd.read_excel('/networks/vtn/lab/SDMC_labscience/assays/MAAR/UW/SDMC_materials/MAAR_info.xlsx')
	maar_metadata = maar_metadata.iloc[:9]
	maar_metadata = maar_metadata.rename(columns={'Name in incoming data':'assay_subtype_lab'})

	metadata_usecols = [
	    'assay_name',
	    'assay_subtype',
	    'data_precision',
	    'instrument',
	    'instrument_serial',
	    'lab_software',
	    'assay_subtype_lab',
	    'LLOD',
	]

	# clean
	maar_metadata.LLOD = maar_metadata.LLOD.str.replace("copies/mL","").str.strip()

	df = df.merge(
	    maar_metadata[metadata_usecols], on='assay_subtype_lab', how='outer'
	)

	df = df.drop(columns='assay_subtype_lab')
	df = df.rename(columns={'guspec1':'guspec'}) # sdmc_tools expects a column called guspec. need to fix.
	df = df.rename(columns={'instrument_serial':'instrument_serialno'})

	ldms = access_ldms.pull_one_protocol('hvtn', 206)
	# ldms = ldms.loc[(ldms.guspec.isin(df.guspec))]

	md = {
	    'network':'HVTN',
	    'specrole':'Sample',
	    'upload_lab_id':'UW',
	    'assay_lab_name':'University of Washington Virology'
	}

	outputs = sdmc_tools.standard_processing(
	    input_data=df,
	    input_data_path=input_path,
	    guspec_col='guspec',
	    network='hvtn',
	    metadata_dict=md,
	    ldms=ldms,
	)

	# sdmc_tools expects a column called guspec. need to fix.
	outputs = outputs.rename(columns={'guspec':'guspec1'})

	# check lab-submitted ptid against ldms, then drop
	assert (outputs.ptid.astype(int)!=outputs.pid.astype(int)).sum()==0
	outputs.ptid = outputs.ptid.astype(int)
	outputs = outputs.drop(columns='pid')

	# switch format from guspec1/2/3 to guspec_core + guspec_aliquots
	outputs_updated = outputs.copy()
	outputs_updated['guspec_core'] = outputs_updated.guspec1.str.rpartition("-")[0]
	outputs_updated['guspec_aliquots'] = outputs_updated.guspec1.str.rpartition("-")[2]+", "+outputs_updated.guspec2.str.rpartition("-")[2]+", "+outputs_updated.guspec3.str.rpartition("-")[2]

	reorder = [
	    'network',
	    'protocol',
	    'guspec_core',
	    'guspec_aliquots',
	    'specrole',
	    'upload_lab_id',
	    'assay_lab_name',
	    'ptid',
	    'visitno',
	    'drawdt',
	    'spectype',
	    'spec_primary',
	    'spec_additive',
	    'spec_derivative',
	    'assay_name',
	    'assay_subtype',
	    'lab_software',
	    'instrument',
	    'instrument_serialno',
	    'result_qualitative',
	    'result_quantitative',
	    'result_detail',
	    'result_units',
	    'llod',
	    'data_precision',
	    'sdmc_processing_datetime',
	    'sdmc_data_receipt_datetime',
	    'input_file_name',
	]

	set(reorder).symmetric_difference(outputs_updated.columns)
	outputs_updated = outputs_updated[reorder]

	# ensure they share the same guspec core
	gus = outputs[[i for i in outputs.columns if 'guspec' in i]]\
	gus.guspec1 = gus.guspec1.str.rpartition("-")[0]
	gus.guspec2 = gus.guspec2.str.rpartition("-")[0]
	gus.guspec3 = gus.guspec3.str.rpartition("-")[0]

	assert (gus.guspec1!=gus.guspec2).sum() == 0
	assert (gus.loc[gus.guspec3.notna()].guspec1!=gus.loc[gus.guspec3.notna()].guspec3).sum()==0

	# outputs.lab_software.unique()

	outputs.llod = outputs.llod.astype(float)

	today = datetime.date.today().isoformat()
	savedir = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN206/assays/MAAR/misc_files/data_processing/'
	outputs.to_csv(savedir + f"HVTN206_MAAR_processed_{today}.txt", sep='\t', index=False)

	# check against manifest from nick
	manifest = pd.read_csv(
	    '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN206/assays/MAAR/misc_files/manifest_512-015-2025000193.txt',
	    sep='\t'
	)

	long_outputs = outputs.melt(
	    id_vars=list(set(outputs.columns).difference(['guspec1','guspec2','guspec3'])),
	    value_vars=['guspec1','guspec2','guspec3'],
	    value_name='guspec',
	    var_name='aliquot'
	)
	long_outputs = long_outputs.loc[long_outputs.guspec.notna()]

	missing_from_data = set(manifest.GLOBAL_ID).difference(long_outputs.guspec)
	missing_from_data

	missing_from_manifest = set(long_outputs.guspec).difference(manifest.GLOBAL_ID)
	missing_from_manifest

	outputs.guspec1.str.rpartition("-")[2].unique()
	outputs.guspec2.str.rpartition("-")[2].unique()
	outputs.guspec3.str.rpartition("-")[2].unique()

	long_outputs.loc[long_outputs.guspec.str.contains("0379-00MFYC00"),['guspec']].drop_duplicates()
	long_outputs['guspec_core'] = long_outputs.guspec.str.rpartition("-")[0]

	cores_missing_from_data = [i.rpartition("-")[0] for i in missing_from_data]
	cores_missing_from_manifest = [i.rpartition("-")[0] for i in missing_from_manifest]

	# guspecs missing from both manifest and data
	set(manifest.GLOBAL_ID.str.rpartition("-")[0]).difference(long_outputs.guspec_core)
	long_outputs.loc[long_outputs.guspec_core.isin(cores_missing_from_data)].groupby(['guspec_core','guspec']).count()[['network']]
	long_outputs.loc[long_outputs.guspec_core.isin(cores_missing_from_manifest)].groupby(['guspec_core','guspec']).count()[['network']]

	# check against manifest from lab
	lab_manifest = pd.read_excel(
	    '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN206/assays/MAAR/misc_files/lab_submission/lab_submitted_manifest/HVTN206-v2.0.xlsx',
	)

	missing_from_data = set(lab_manifest['Global Spec ID']).difference(long_outputs.guspec)

	lab_manifest['guspec_core'] = lab_manifest['Global Spec ID'].str.rpartition("-")[0]

	# these are the samples that weren't run for RNA
	not_dones = outputs.loc[(outputs.assay_subtype=='Abbott HIV-1 RNA') & (outputs.result_quantitative=="Not Done")].guspec1.unique().tolist()

	# confirm samples in manifest missing from data are exactly the expected ones
	assert set([i.rpartition("-")[0] for i in not_dones]).symmetric_difference([i.rpartition("-")[0] for i in missing_from_data])==set()

	melt_cols = [
	    'network',
	    'protocol',
	    'ptid',
	    'visitno',
	    'drawdt',
	    'spectype',
	    'spec_primary',
	    'spec_additive',
	    'spec_derivative',
	    'assay_name',
	    'assay_subtype',
	    'lab_software',
	    'instrument',
	    'instrument_serialno',
	    'result_qualitative',
	    'result_quantitative',
	    'result_detail',
	    'result_units',
	    'llod',
	    'data_precision',
	    'sdmc_processing_datetime',
	    'sdmc_data_receipt_datetime',
	    'input_file_name',
	]

	# pivot summaries
	summary = outputs.melt(
	    id_vars=melt_cols,
	    value_vars=['guspec1','guspec2','guspec3'],
	    value_name='guspec'
	)

	summary = summary.loc[summary.guspec.notna()]

	summary1 = pd.pivot_table(
	    summary,
	    index=['visitno','ptid'],
	    columns='assay_subtype',
	    aggfunc='count'
	)[['guspec']]
	summary1.to_excel(savedir + "HVTN206_MAAR_sample_summary.xlsx")

	summary2 = outputs.copy()
	summary2['is_reactive'] = False
	summary2.loc[summary2.result_qualitative=="Ab Reactive", 'is_reactive'] = True

	summary2 = summary2.groupby(['visitno','ptid'])[['is_reactive']].sum().reset_index()
	summary2 = summary2.rename(columns={'is_reactive':'count_of_reactives'})
	summary2 = summary2.set_index(['visitno','ptid'])

	summary2.to_excel(savedir + "HVTN206_MAAR_sample_reactivity_summary.xlsx")

if __name__=="__main__":
	main()