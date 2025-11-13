## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 11/12/2025
# Purpose: Blocking
## ---------------------------------------------------------------------------##
import pandas as pd
import numpy as np
import os
import datetime as datetime

import sdmc_tools.process as sdmc_tools
import sdmc_tools.access_ldms as access_ldms

def main():
	# read in blocking data + both columns (columns in two separate rows)
	blocking_path = '/trials/vaccine/p115/s001/qdata/elisa_blocking_pass-through/251001 HVTN 115 and 135 Post 5 Blocking Summary.xlsx'
	blocking_cols = pd.read_excel(blocking_path, sheet_name='Blocking Summary', skiprows=6, nrows=1)
	blocking = pd.read_excel(blocking_path, sheet_name='Blocking Summary', skiprows=7)

	# merge both rows into columns
	bcols = list(blocking.columns)
	bcols0, bcols1 = bcols[:10], bcols[10:]
	bcols11 = [x for x in list(blocking_cols.columns)[10:][::2] for _ in range(2)]
	bcols1 = [f"{j}|{i}" for (i,j) in zip(bcols1,bcols11)]
	blocking.columns = bcols0 + bcols1

	# cast to long on analyte
	blocking = blocking.melt(
	    id_vars=bcols0,
	    value_vars=bcols1,
	    value_name='result',
	    var_name='analyte|metric',
	)

	# split out antigen and metric
	antigen_metrics = blocking['analyte|metric'].str.split("|", expand=True).rename(columns={0:'analyte',1:'metric'})
	blocking = pd.concat([blocking, antigen_metrics], axis=1).drop(columns='analyte|metric')
	blocking.metric = blocking.metric.str.split(".", expand=True)[0]

	# back to wide on metric
	blocking = blocking.pivot(
	    columns='metric',
	    index=['Original ID', 'Date Drawn', 'Group/Treatment', 'Tx', 'Part',
	       'Study ID', 'Pub ID', 'Subject ID', 'Visit', 'Day', 'analyte'],
	    values='result'
	).reset_index()

	# read in control data, update column names
	controls = pd.read_excel(blocking_path, sheet_name='Controls', skiprows=19)
	ccols0 = list(controls.columns)
	ccols1 = list(controls.iloc[0])
	ccols0 = ccols0[:3] + [x for x in ccols0[3:][::2] for _ in range(2)]

	controls_cols = [f"{i}|{j}" for (i,j) in zip(ccols0,ccols1)]
	controls_cols[9:11] = ['CH106 cold x CH103uca|OD', 'CH106 cold x CH103uca|% Blocking']
	controls.columns = controls_cols

	controls = controls.iloc[1:].reset_index(drop=True)
	controls = controls.dropna(how='all', axis=1)

	# cast to long along metric
	controls = controls.melt(
	    id_vars=[
	        'Unnamed: 0|nan',
	        'Standards|mcg/ml',
	        '1|Wells'
	    ],
	    value_vars=[
	        'CH106|OD',
	        'CH106|% Blocking',
	        'CH103uca|OD',
	        'CH103uca|% Blocking',
	        'CH106 cold x CH103uca|OD',
	        'CH106 cold x CH103uca|% Blocking',
	        'sCD4|OD',
	        'sCD4|% Blocking',
	        'CH235.12|OD',
	        'CH235.12|% Blocking',
	        'CH235uca|OD',
	        'CH235uca|% Blocking',
	        'CH01|OD',
	        'CH01|% Blocking',
	        'PGT145|OD',
	        'PGT145|% Blocking',
	    ],
	    var_name='analyte|metric',
	    value_name='result'
	)

	# split out analyte and metric into separate columns
	control_metrics = controls['analyte|metric'].str.split("|", expand=True).rename(columns={0:'analyte',1:'metric'})
	controls = pd.concat([controls, control_metrics], axis=1)
	controls = controls.drop(columns='analyte|metric')

	# cast back to wide along metric
	controls = controls.pivot(
	    values='result',
	    columns='metric',
	    index=['Unnamed: 0|nan', 'Standards|mcg/ml', '1|Wells','analyte']
	).reset_index()

	# read in averaged control data (block at top of tab)
	control_avg = pd.read_excel(blocking_path, sheet_name='Controls', skiprows=4)
	control_avg = control_avg.iloc[1:10].dropna(how='all', axis=1)

	control_avg_cols = [
	    'Standards|mcg/ml cold antibody',
	    'CH106|OD',
	    'CH106|pct_blocking',
	    'CH65 x CH106|OD',
	    'CH65 x CH106|pct_blocking',
	    'CH103uca|OD',
	    'CH103uca|pct_blocking',
	    'CH106 cold x CH103uca|OD',
	    'CH106 cold x CH103uca|pct_blocking',
	    'CH65 x CH103_UCA|OD',
	    'CH65 x CH103_UCA|pct_blocking',
	    'sCD4|OD',
	    'sCD4|pct_blocking',
	    'CH235.12|OD',
	    'CH235.12|pct_blocking',
	    'CH235uca|OD',
	    'CH235uca|pct_blocking',
	    'CH01|OD',
	    'CH01|pct_blocking',
	    'PGT145|OD',
	    'PGT145|pct_blocking',
	]
	control_avg.columns = control_avg_cols

	# cast to long along analyte
	control_avg = control_avg.melt(
	    id_vars=['Standards|mcg/ml cold antibody'],
	    value_vars=list(set(control_avg_cols).difference(['Standards|mcg/ml cold antibody'])),
	    value_name='result',
	    var_name='analyte|metric'
	)

	# split analyte and metric into separate columns 
	cavg_antigen_metrics = control_avg['analyte|metric'].str.split("|", expand=True).rename(columns={0:'analyte',1:'metric'})
	control_avg = pd.concat([control_avg.drop(columns='analyte|metric'), cavg_antigen_metrics], axis=1)

	# subset to relevant column, cast back to wide on metric
	control_avg = control_avg[['Standards|mcg/ml cold antibody',
	       'result', 'analyte', 'metric']]
	control_avg = control_avg.pivot(
	    values='result',
	    columns='metric',
	    index=['analyte','Standards|mcg/ml cold antibody'],
	).reset_index()

	# rename columns and merge
	control_avg.columns = ['analyte','control_concentration','result_average','result_pct_blocking_avg']
	controls.columns = ['replicate','control_concentration','well','analyte','result_pct_blocking','result']
	controls.replicate = controls.replicate.astype(int)
	controls = controls.merge(control_avg, on=['analyte','control_concentration'], how='outer').sort_values(by=['analyte','control_concentration','replicate'])


	# drop rows with no results (these were added by the outer merge)
	controls = controls.loc[(
	    (controls.result_pct_blocking.notna()) |
	    (controls.result.notna()) |
	    (controls.result_average.notna()) |
	    (controls.result_pct_blocking_avg.notna())
	)]

	# standardize analyte names
	antigen_map = {
	    'CH01': 'CH01',
	    'CH103uca': 'CH103_UCA',
	    'CH106': 'CH106',
	    'CH106 cold x CH103uca': 'CH106 cold x CH103uca',
	    'CH235.12': 'CH235.12',
	    'CH235uca': 'CH235_UCA',
	    'CH65 x CH103_UCA': 'CH65 x CH103_UCA',
	    'CH65 x CH106': 'CH65 x CH106',
	    'PGT145': 'PGT145',
	    'sCD4': 'sCD4',
	}
	controls.analyte = controls.analyte.map(antigen_map)

	antigen_vals = pd.DataFrame({'Controls':antigen_map.keys(), 'Control Averages':antigen_map.values()})
	antigen_vals = antigen_vals.loc[~antigen_vals.Controls.str.contains("x")]
	antigen_vals['Samples'] = ['CH01','CH103_UCA','CH106','CH235.12*','CH235_UCA*','PGT145','sCD4']

	assert set(blocking.analyte).symmetric_difference(antigen_vals.Samples) == set()

	blocking.analyte = blocking.analyte.map(antigen_vals.set_index('Samples')['Control Averages'])

	# rename column
	blocking = blocking.rename(columns={'% Block':'result_pct_blocking', 'Stdev':'result_pct_blocking_stddev'})

	# add specrole and concat
	controls['specrole'] = 'Standard'
	blocking['specrole'] = 'Sample'
	blocking = pd.concat([blocking, controls])

	assert blocking['Original ID'].isna().sum() == len(controls)

	blocking = blocking.rename(columns={'Original ID':'guspec'})

	# Assays conducted on CH505TFchim.6R.SOSIP.664v4.1_10lnQQ-avi-BIO/293F 
	# except "*" conducted on CH505M5.G458Ychim.6R.SOSIP.MD39_2P_csorta.10lnQQ-avi-Bio/GnTI-
	ag_map = {i:"CH505TFchim.6R.SOSIP.664v4.1_10lnQQ-avi-BIO/293F" for i in [j for j in blocking.analyte.unique() if 'x' not in j]}
	ag_map['CH235.12'] = 'CH505M5.G458Ychim.6R.SOSIP.MD39_2P_csorta.10lnQQ-avi-Bio/GnTI-'
	ag_map['CH235uca'] = 'CH505M5.G458Ychim.6R.SOSIP.MD39_2P_csorta.10lnQQ-avi-Bio/GnTI-'

	blocking['target'] = blocking.analyte.map(ag_map)

	ldms = access_ldms.pull_multiple_protocols('hvtn', [115, 135])
	ldms = ldms.loc[ldms.guspec.isin(blocking.guspec.tolist())]

	# standard processing
	md = {
	    'upload_lab_id': 'BH',
	    'assay_lab_name':'Haynes Lab (Duke)',
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

	# ensure these match
	outputs[['visit','visitno']].drop_duplicates()

	# visual check -- these should match
	outputs.loc[outputs.ptid.astype(float)!=outputs.subject_id.astype(float),['ptid','subject_id']].drop_duplicates()

	# visual check -- these should match
	outputs.loc[outputs.drawdt!=outputs.date_drawn,['drawdt','date_drawn']].drop_duplicates()

	# visual check -- these should match
	outputs[['study_id','protocol']].drop_duplicates()

	# everything matched, can drop lab-submitted metadata
	outputs = outputs.drop(columns=['subject_id','visit','date_drawn','study_id'])

	# lab shouldn't be the source for this type of info
	outputs = outputs.drop(columns=[
	    'pub_id',
	    'day',
	    'tx',
	    'part',
	    'group/treatment',
	])


	# network should only be associated with sample rows
	outputs.loc[outputs.specrole=='Sample','network'] = 'HVTN'

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
	    'control_concentration',
	    'control_concentration_units',
	    'result',
	    'result_average',
	    'result_units',
	    'result_pct_blocking',
	    'result_pct_blocking_avg',
	    'result_pct_blocking_stddev',
	    'well',
	    'replicate',
	    'sdmc_processing_datetime',
	    'sdmc_data_receipt_datetime',
	    'input_file_name'
	]

	assert set(reorder).symmetric_difference(outputs.columns) == set()
	outputs = outputs[reorder]

	# cast to numeric
	outputs.protocol = outputs.protocol.astype(float)

	# visual check
	outputs[['protocol','specrole']].drop_duplicates()

	# visual check
	outputs.loc[outputs.guspec.isna(),['protocol','specrole']].drop_duplicates()

	# split out by protocol + controls
	outputs115 = outputs.loc[(outputs.protocol==115.) | (outputs.specrole=='Standard')]
	outputs135 = outputs.loc[(outputs.protocol==135.) | (outputs.specrole=='Standard')]

	# save by protocol
	today = datetime.date.today().isoformat()

	# savedir115 = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN115/assays/blocking_ELISA/misc_files/data_processing/'
	# outputs115.to_csv(savedir115 + f'HVTN115_ELISA_bnab_blocking_processed_{today}.txt', sep="\t", index=False)

	# savedir135 = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN135/assays/blocking_ELISA/misc_files/data_processing/'
	# outputs135.to_csv(savedir135 + f'HVTN135_ELISA_bnab_blocking_processed_{today}.txt', sep="\t", index=False)

	# # save pivot summaries
	# summary115 = pd.pivot_table(
	#     outputs115,
	#     index=['specrole','visitno','ptid'],
	#     columns=['analyte'],
	#     values='sdmc_processing_datetime',
	#     aggfunc='count',
	#     dropna=False
	# )
	# summary115.dropna(how='all').to_excel(savedir115 + "HVTN115_ELISA_bnab_blocking_summary.xlsx")

	# summary135 = pd.pivot_table(
	#     outputs135,
	#     index=['specrole','visitno','ptid'],
	#     columns=['analyte'],
	#     values='sdmc_processing_datetime',
	#     aggfunc='count',
	#     dropna=False
	# )
	# summary135.dropna(how='all').to_excel(savedir135 + "HVTN135_ELISA_bnab_blocking_summary.xlsx")

	# check against manifest
	manifest1 = pd.read_csv(
	    '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN115/assays/binding_ELISA/misc_files/manifests/481-373-0000043537.txt',
	    sep='\t'
	)
	manifest2 = pd.read_csv(
	    '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN115/assays/binding_ELISA/misc_files/manifests/512--2025000186.txt',
	    sep='\t'
	)

	manifest = pd.concat([manifest1, manifest2])

	# looks like some of them are in the manifest
	set(outputs135.guspec).intersection(manifest.GLOBAL_ID)

	# but a lot of them arent
	set(outputs135.guspec).symmetric_difference(manifest.GLOBAL_ID)

	# 115 looks good
	set(outputs115.guspec).symmetric_difference(manifest.loc[manifest.PROTOCOL==115.].GLOBAL_ID)

if __name__=="__main__":
	main()