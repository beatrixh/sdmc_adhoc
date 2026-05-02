## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 04/30/2026
# Purpose:  run ad hoc processing on ju's retagged air files
# Note: 	unlike usual, we arent touching ju's columns at all.
#			only merging on ldms + metadata and running ldms checks
## ---------------------------------------------------------------------------##
import pandas as pd
import numpy as np
import sdmc_tools.access_ldms as access_ldms
import sdmc_tools.process as sdmc
import datetime
import os

## process long
def process_dataset(input_data_path):
	df = pd.read_csv(
	    input_data_path,
	    sep="\t"
	)
	df.Global_Spec_ID = df.Global_Spec_ID.str.replace("\t","")

	# rerun ldms checks
	ldms = access_ldms.pull_one_protocol('hvtn', 807)
	ldms['guspec_core'] = ldms.guspec.str.rpartition("-")[0]
	ldms = ldms.drop(columns='guspec').drop_duplicates()

	df_check = df[['Global_Spec_ID','Visit','PTID']].drop_duplicates().rename(columns={'PTID':'ptid_lab','Global_Spec_ID':'guspec_core','Visit':'visitno_lab'})
	df_check.guspec_core = df_check.guspec_core.str.replace("\t","")

	df_check = df_check.merge(ldms, on='guspec_core', how='left')
	df_check = df_check.rename(columns={'txtpid':'ptid_ldms', 'vidval':'visitno_ldms'})

	assert len(df_check.loc[(df_check.visitno_lab.astype(float) != df_check.visitno_ldms.astype(float)),['guspec_core','visitno_lab','visitno_ldms']]) == 0
	assert len(df_check.loc[(df_check.ptid_lab.astype(float)!=df_check.ptid_ldms.astype(float)),['guspec_core','ptid_lab','ptid_ldms']]) == 0

	# make this guspec for merge code
	ldms = ldms.rename(columns={'guspec_core':'guspec'})

	md = {
	    'network':'HVTN',
	    'specrole':'Sample',
	    'upload_lab_id':'O1',
	    'assay_lab_name':'McElrath Lab',
	    'assay_type':'BCR Sequencing',
	    'assay_subtype':'10x',
	    'assay_details':'10x Chromium GEM-X Single Cell 5’ v3 Sequencing',
	    'assay_precision':'Quantitative',
	    'instrument':'NovaSeq X+ (Illumina)',
	    'lab_software':'proSCessoR, Cell Ranger, V(D)J Reference OGRDB',
	    'lab_software_version':'proSCessoR v0.3.0, Cell Ranger v10.0.0, V(D)J Reference OGRDB 2026.01.26',
	    'sdmc_bioinformatics_software':'BALDR',
	    'sdmc_bioinformatics_software_version':'v0.7.8 (commit hash: 7febff8)',
	}

	outputs = sdmc.standard_processing(
	    df,
	    input_data_path=input_data_path,
	    guspec_col='Global_Spec_ID',
	    network='HVTN',
	    metadata_dict=md,
	    ldms=ldms,
	    cols_to_lower=False,
	)
	outputs = outputs.rename(columns={'guspec':'guspec_core'})

	first_vars = [
	    'network',
	    'protocol',
	    'guspec_core',
	    'specrole',
	    'visitno',
	    'ptid',
	    'drawdt',
	    'spectype',
	    'spec_primary',
	    'spec_additive',
	    'spec_derivative',
	    'upload_lab_id',
	    'assay_lab_name',
	    'assay_type',
	    'assay_subtype',
	    'assay_details',
	    'assay_precision',
	    'instrument',
	    'lab_software',
	    'lab_software_version',
	    'sdmc_bioinformatics_software',
	    'sdmc_bioinformatics_software_version',
	    'input_file_name',
	    'sdmc_data_receipt_datetime',
	    'sdmc_processing_datetime',
	]

	reorder = first_vars + df.columns.tolist()
	assert set(reorder).symmetric_difference(outputs.columns) == set()

	outputs = outputs[reorder]
	outputs = outputs.drop(columns=['Visit','PTID','Global_Spec_ID'])

	outputs.columns = [i.lower().replace(" ","_") for i in outputs.columns]

	return outputs

today = datetime.date.today().isoformat()
savedir = "/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN807/assays/BCR_sequencing/misc_files/data_processing/"


input_data_path="/trials/vaccine/p807/s001/qdata/LabData/BCR_sequencing_pass-through/uploaded_by_lab/Hyrien/20260428-02/20260427_hvtn807_airr.tsv"
outputs_long = process_dataset(input_data_path)
outputs_long.to_csv(
    savedir + f"HVTN807_BCR_sequencing_regenerated_AIRR_long_processed_{today}.txt",
    sep="\t", index=False
)

input_data_path="/trials/vaccine/p807/s001/qdata/LabData/BCR_sequencing_pass-through/uploaded_by_lab/Hyrien/20260428-02/20260427_hvtn807_wide.tsv"
outputs_wide = process_dataset(input_data_path)
outputs_wide.to_csv(
    savedir + f"HVTN807_BCR_sequencing_regenerated_AIRR_wide_processed_{today}.txt",
    sep="\t", index=False
)