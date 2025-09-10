## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 09/08/2025
# Purpose: Ad hoc processing of covpn5001 proteomics data
## ---------------------------------------------------------------------------##
import pandas as pd
import numpy as np
import os
import datetime as datetime

import sdmc_tools.process as sdmc_tools
import sdmc_tools.access_ldms as access_ldms

def main():
    # read in data
    input_data_path = '/trials/covpn/p5001/s001/qdata/LabData/proteomics_pass-through/Postcleanup_COVPN_Olink_NPX_Matrix_merged_metadata.csv'
    data = pd.read_csv(input_data_path)
    data = data.drop(columns = 'Unnamed: 0')

    # cast to long
    analytes = data.columns[1:-60]
    demos = data.columns.tolist()[-60:]
    cols = list(data.columns)
    id_vars = ['UniqueSampleID'] + demos

    data = pd.melt(
        data,
        id_vars=id_vars,
        value_vars=analytes,
        var_name='analyte',
        value_name='result'
    )

    data = data.rename(columns={'external.subject.ID':'ptid'})

    # read in manifest to get ptid/guspec map
    misc_files_dir = '/networks/vtn/lab/SDMC_labscience/studies/CoVPN/CoVPN5001/assays/proteomics/misc_files/'
    manifest1 = pd.read_csv(misc_files_dir + '481-999110-0000051745.txt', sep="\t")
    manifest2 = pd.read_csv(misc_files_dir + '481-999110-0000051746.txt', sep="\t")
    manifest3 = pd.read_csv(misc_files_dir + '512--2023000097.txt', sep="\t")

    manifest = pd.concat([manifest1, manifest3])
    manifest = manifest.loc[manifest.DER=='SER']
    manifest = manifest.rename(columns={'PID':'ptid', 'GLOBAL_ID':'guspec'})

    # confirm one guspec per sample
    assert len(manifest) == manifest.ptid.nunique()

    data = data.merge(manifest[['ptid','guspec','SHIPPED_FROM']], on='ptid', how='left')

    # merge on metadata from barcode mapping file they shared
    barcode_mapping_path = '/networks/vtn/lab/SDMC_labscience/studies/CoVPN/CoVPN5001/assays/proteomics/lab_submitted_data/CoVPN5001 proteomic internal_external barcode mapping.xlsx'
    barcode = pd.read_excel(barcode_mapping_path, sheet_name='data', skiprows=1)

    # make sure these are a match
    assert len(set(data.UniqueSampleID).symmetric_difference(barcode['AIFI Barcode'])) == 0
    assert barcode['AIFI Barcode'].nunique() == len(barcode)

    data = data.merge(barcode, left_on='UniqueSampleID', right_on='AIFI Barcode', how='outer')
    data = data.rename(columns={
        'Type (cntp_name)': 'spectype_from_barcode',
        '<SLIMSGUID>': 'slims_gu_id',
        'AIFI Barcode': 'lab_internal_sample_id',
        'Subject': 'lab_internal_subject_id',
        'Cohort': 'lab_internal_cohort_id',
    })

    # some of the cols we've merged on are now redundant
    assert (data.ptid!=data['External Subject ID']).sum() == 0
    assert (data['internal.subject']!=data.lab_internal_subject_id).sum() == 0

    dropcols = [
        'UniqueSampleID',
        'internal.subject',
        'External Subject ID',
    ]
    data = data.drop(columns=dropcols)

    # standard processing
    ldms = access_ldms.pull_one_protocol('covpn', 5001)
    ldms = ldms.loc[ldms.guspec.isin(data.guspec)]

    md = {
        'network':'CoVPN',
        'upload_lab_id': 'N/A',
        'assay_lab_name':'Allen Institute',
        'assay_type': 'Proteomics',
        'specrole': 'Sample',
        'instrument':'Olink Explore',
        'lab_software':'NPX Map',
        'result_units': 'NPX',
    }

    outputs = sdmc_tools.standard_processing(
        input_data=data,
        input_data_path=input_data_path,
        guspec_col='guspec',
        network='CoVPN',
        metadata_dict=md,
        ldms=ldms,
    )


    # drop this (at least for now; unless the lab makes it sound useful)
    assert (outputs.spectype_from_barcode!=outputs.spectype).sum() == 0
    outputs[['visit','visitno']].drop_duplicates()

    assert (outputs.ptid_x.astype(int)!=outputs.ptid_y.astype(int)).sum() == 0
    outputs = outputs.drop(columns=[
        'visit',
        'spectype_from_barcode',
        'ptid_y',
        'studyid',
    ])
    outputs = outputs.rename(columns={
        'ptid_x':'ptid'
    })

    # reorder columns
    reorder = [
        'network',
        'protocol',
        'specrole',
        'guspec',
        'pubid',
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
        'instrument',
        'lab_software',
        'analyte',
        'result',
        'result_units',
        'site',
        'country',
        'region',
        'age',
        'sex',
        'race',
        'racesp',
        'ethnic',
        'weightbl',
        'heightbl',
        'bmibl',
        'syfl',
        'symstdy',
        'covstdy',
        'covposdy',
        'termreas',
        'termfl',
        'enrgr1',
        'enrgr1n',
        'maxvisit',
        'v6fl',
        'v7fl',
        'mhhiv',
        'mhhyper',
        'mhchf',
        'mhcea',
        'mhdiab',
        'mhdiabdx',
        'mhckd',
        'mhisd',
        'mhany',
        'smokenow',
        'smokever',
        'pregfl',
        'ppfl',
        'cohort',
        'mhcov',
        'cortfl',
        'dexamfl',
        'hospadm',
        'icu',
        'oxy_any',
        'm2_oxy',
        'mxsevsr',
        'mxsev',
        'mxsevox',
        'm2_syfl',
        'm2_syfl_clean',
        'first_clearance_cdays',
        'first_clearance_category',
        'vaccfl',
        'vacb4onst',
        'vacstatus',
        'vaccseries',
        'lineage',
        'who_label',
        'lab_internal_cohort_id',
        'lab_internal_sample_id',
        'lab_internal_subject_id',
        'slims_gu_id',
        'shipped_from',
        'sdmc_processing_datetime',
        'sdmc_data_receipt_datetime',
        'input_file_name'
    ]
    assert len(set(reorder).symmetric_difference(outputs.columns)) == 0
    outputs = outputs[reorder]

    # save outputs
    savedir = "/networks/vtn/lab/SDMC_labscience/studies/CoVPN/CoVPN5001/assays/proteomics/misc_files/data_processing/"
    today = datetime.date.today().isoformat()
    outputs.to_csv(savedir + f"CoVPN5001_Allen_Institute_Proteomics_Processed_{today}.txt", index=False, sep="\t")

    # qdata_dir = input_data_path.rpartition("/")[0]
    # os.listdir(qdata_dir)

    # pivot summary
    summary = pd.pivot_table(
        outputs,
        index='analyte',
        columns='ptid',
        aggfunc='count',
        fill_value=0
    )[['result']]
    summary.to_excel(savedir + "covpn5001_allen_inst_proteomics_pivot_summary.xlsx")

    summary2 = pd.pivot_table(
        outputs,
        index='ptid',
        columns='visitno',
        aggfunc='count',
        fill_value=0
    )[['result']]

    summary2.to_excel(savedir + "covpn5001_allen_inst_proteomics_ptid_visitno_pivot_summary.xlsx")
if __name__=="__main__":
    main()
