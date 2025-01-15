## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 08/26/2024
# Purpose:  Ad hoc data processing of PVMA data from Fouda Lab
## ---------------------------------------------------------------------------##
import pandas as pd
import numpy as np
import datetime
import yaml
import os

import sdmc_tools.process as sdmc
import sdmc_tools.constants as constants
## ---------------------------------------------------------------------------##

def main():
    input_data_path = '/trials/vaccine/p135/s001/qdata/LabData/PVMA_pass-through/20240822_HVTN135_DataSummary_V2.xlsx'
    data = pd.read_excel(input_data_path)

    rename = {
        'GLOBAL_ID':'guspec',
        'Tube Label (WCM Internal Only)': 'lab_internal_sample_id',
        'Plate #': 'plate_number',
        'internal_qc':'lab_internal_qc'
    }
    data = data.rename(columns=rename)

    result_cols = [i for i in data.columns if '(' in i]
    non_result_cols = list(set(data.columns).difference(result_cols))
    data = data.melt(id_vars = non_result_cols, value_vars = result_cols, var_name='antigen_string', value_name='result')

    results = data.antigen_string.str.rpartition(" ", expand=True)
    results.columns = ['antigen', 'na', 'beadset_num']
    data = pd.concat([data, results], axis=1)
    data.beadset_num = data.beadset_num.str[1:-1].astype(int)

    data.columns = [i.lower().replace(" ","_") for i in data.columns]

    ldms = pd.read_csv(
        constants.LDMS_PATH_HVTN,
        usecols=constants.STANDARD_COLS,
        dtype=constants.LDMS_DTYPE_MAP
    )
    ldms = ldms.loc[ldms.lstudy==135.]
    ldms = ldms.loc[ldms.guspec.isin(data.guspec)]

    md = {
        'network': 'HVTN',
        'upload_lab_id': 'P5',
        'assay_lab_name': 'Fouda Lab (Cornell)',
        'instrument': 'Bio-Rad BioPlex 200',
        'lab_software_version': 'BioPlex 6.2',
        'assay_type': 'BAMA',
        'assay_subtype': 'PVMA',
        'specrole': 'Sample',
        'result_units': 'RLU'
    }

    outputs = sdmc.standard_processing(
        input_data=data,
        input_data_path=input_data_path,
        guspec_col="guspec",
        network="HVTN",
        metadata_dict=md,
        ldms=ldms
    )

    drop_cols = ['assay_lab', 'comments', 'der', 'upload_date', 'antigen_string', 'na']
    outputs = outputs.drop(columns=drop_cols)

    # rename to match the ELISA'd HepB data
    outputs = outputs.rename(columns={'sample_dilution':'dilution'})

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
        'instrument',
        'lab_software_version',
        'antigen',
        'beadset_num',
        'dilution',
        'result',
        'result_units',
        'assay_date',
        'assay_tech',
        'internal_qc',
        'lab_internal_sample_id',
        'plate_number',
        'sdmc_processing_datetime',
        'sdmc_data_receipt_datetime',
        'input_file_name'
    ]
    mismatch = set(outputs.columns).symmetric_difference(reorder)
    if len(mismatch) > 0:
        raise Exception("Trying to add or drop columns!")
    outputs = outputs[reorder]

    outputs.assay_date = outputs.assay_date.astype(str)
    # i think excel tried to make this a time
    outputs['dilution'] = outputs['dilution'].astype(str).str[1:5]

    # drop hepB; they're going to assay it separately
    outputs = outputs.loc[outputs.antigen!="HepB"]

    savedir = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN135/assays/bama_pvma/misc_files/data_processing/'
    today = datetime.date.today().isoformat()
    outputs.to_csv(savedir + f"DRAFT_HVTN135_PVMA_FOUDA_processed_{today}.txt", sep="\t", index=False)

    ## generate pivot summary
    pivot_summary = pd.pivot_table(data=outputs, index='ptid', columns='visitno', aggfunc='count')['result']
    pivot_summary.to_excel(savedir + "HVTN135_PVMA_ptid_visitno_summary.xlsx")

def checks():
    savedir = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN135/assays/bama_pvma/misc_files/data_processing/'
    files = os.listdir(savedir)
    fname = [i for i in files if 'HVTN135_PVMA_FOUDA_processed' in i][0]
    output = pd.read_csv(savedir + fname, sep="\t")

    manifest_path = "/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN135/assays/bama_pvma/misc_files/HVTN 135 v12 Shipment to Cornell_Fouda.xlsx"
    manifest = pd.read_excel(manifest_path, sheet_name="Sheet2")
    manifest_mismatch = set(outputs.guspec).symmetric_difference(manifest.GLOBAL_ID)
    if len(manifest_mismatch) > 0:
        raise Exception("Manifest doesn't match data received")

    original_data_path = "/trials/vaccine/p135/s001/qdata/LabData/PVMA_pass-through/preliminary_data/April2024_HVTN 135 Summary Data_V1.xlsx"
    og_data = pd.read_excel(original_data_path)
    og_mismatch = set(outputs.guspec).symmetric_difference(og_data.GLOBAL_ID)
    if len(og_mismatch) > 0:
        raise Exception("Samples don't match prelimiary data")

if __name__=="__main__":
    main()
    checks()
