## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 11/19/2025
# Purpose: Ad hoc processing of Accelevir CA-RNA data for HVTN805
## ---------------------------------------------------------------------------##
## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 11/19/2025
# Purpose: Accelevir
## ---------------------------------------------------------------------------##
import pandas as pd
import numpy as np
import os
import datetime as datetime

import sdmc_tools.process as sdmc_tools
import sdmc_tools.access_ldms as access_ldms


def main():
	carna_path = '/trials/vaccine/p805/s001/qdata/LabData/CA-RNA_pass-through/20251113 - ACDX CARNAv1 Data Report - HVTN805 cumulative.xlsx'
	carna = pd.read_excel(carna_path)

	carna.columns = [i.lower().replace(" ","_").replace("_(carna)","") for i in carna.columns]

	guspecs = carna.global_id.str.split(":", expand=True)
	guspecs.columns = ['guspec','guspec2','guspec3']
	carna = pd.concat([guspecs, carna], axis=1)

	col_rename = {
	    'global_id': 'global_id',
	    'pbmc_count_at_thaw\xa0': 'cell_count_at_thaw',
	    'pbmc_viability': 'cell_viability',
	    'cells_analyzed_per_assay': 'cells_analyzed',
	    'pol_ca-rna_copies/1e6_pbmcs': 'result_pol',
	    'polya_ca-rna_copies/1e6_pbmcs': 'result_polya', 
	    'sample_comments':'comments'   
	}
	carna = carna.rename(columns=col_rename)

	# boolean remap
	carna.sample_quality_issue = carna.sample_quality_issue.map({"N":False})
	carna.useit = carna.useit.map({"Y":True})
	carna.loc[carna.deviation_of_testing.isna(),'deviation_of_testing'] = False

	ldms = access_ldms.pull_one_protocol('hvtn', 805)
	assert set(carna.guspec).union(carna.guspec2).union(carna.guspec3).difference(ldms.guspec).difference({None}) == set()
	guspecs_all = list(set(carna.guspec).union(carna.guspec2).union(carna.guspec3).difference({None}))
	ldms = ldms.loc[ldms.guspec.isin(guspecs_all)]

	metadata = pd.read_excel('/networks/vtn/lab/SDMC_labscience/assays/CA-RNA/Accelevir/SDMC_materials/CA-RNA_info.xlsx', header=None)

	# we'll populate these ourselves
	carna = carna.drop(columns=[
	    'assay_lab_name',
	    'assay_name',
	])

	# standard processing
	md = metadata.set_index(0)[1].to_dict()
	md['Specrole'] = 'Sample'

	outputs = sdmc_tools.standard_processing(
	    input_data=carna,
	    input_data_path=carna_path,
	    guspec_col='guspec',
	    network='hvtn',
	    metadata_dict=md,
	    ldms=ldms
	)

	outputs.collection_date_time = pd.to_datetime(outputs.collection_date_time).dt.date.astype(str)
	outputs.drawdt = outputs.drawdt.astype(str)

	outputs.ptid = outputs.ptid.astype(int)
	outputs.visitno = outputs.visitno.astype(float)
	outputs.protocol = outputs.protocol.astype(float)

	assert (outputs.ptid!=outputs.patient_id).sum() == 0
	assert (outputs.visitno!=outputs.visit_id).sum() == 0
	assert (outputs.collection_date_time!=outputs.drawdt).sum() == 0

	assert len(outputs[['protocol','study_id']].drop_duplicates()) == 1

	outputs = outputs.drop(columns=[
	    'study_id',
	    'patient_id',
	    'visit_id',
	    'collection_date_time',
	    'spec_type',
	    'global_id',
	])

	outputs = outputs.rename(columns={'guspec':'guspec1'})

	reorder = [
	    'network',
	    'protocol',
	    'specrole',
	    'guspec1',
	    'guspec2',
	    'guspec3',
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
	    'assay_date',
	    'cell_count_at_thaw',
	    'cell_viability',
	    'cells_analyzed',
	    'result_pol',
	    'result_polya',
	    'result_units',
	    'comments',
	    'sample_quality_issue',
	    'deviation_of_testing',
	    'useit',
	    'sdmc_processing_datetime',
	    'sdmc_data_receipt_datetime',
	    'input_file_name',
	]

	print(set(outputs.columns).symmetric_difference(reorder))
	assert set(outputs.columns).symmetric_difference(reorder) == set()
	outputs = outputs[reorder]

	savedir = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN805_HPTN093/assays/CA-RNA/misc_files/data_processing/'
	today = datetime.date.today().isoformat()
	outputs.to_csv(
	    savedir + f"HVTN805_Accelevir_CA-RNA_Processed_{today}.txt",
	    sep="\t",
	    index=False
	)

	summary = pd.pivot_table(
	    outputs,
	    index=['ptid'],
	    columns='visitno',
	    aggfunc='count',
	    values='result_pol',
	    fill_value=0
	)

	summary.to_excel(
	    savedir + "HVTN805_Accelevir_CA-RNA_summary.xlsx"
	)

	# check against manifest
	manifest = pd.read_excel(
	    '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN805_HPTN093/assays/IPDA/misc_files/MANIFEST_SHIP_REQ_0352.xlsx',
	)

	set(manifest.GLOBAL_ID).symmetric_difference(guspecs_all)

if __name__=="__main__":
	main()