## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 01/14/2025
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
    input_data_path = '/trials/vaccine/p135/s001/qdata/LabData/HepB_ELISA_pass-through/HVTN_HepB ELISA Binding_Summary_Upload_V2.xlsx'
    df = pd.read_excel(input_data_path)

    rename = {
        'SampleID (Fouda Lab Only)': 'lab_internal_sample_id',
        'GLOBAL_ID': 'guspec',
        'DER': 'DER',
        'Assay Lab': 'Assay Lab',
        'Assay Date': 'assay_date',
        'Plate #': 'plate_no',
        'Corresponding Dilution': 'dilution',
        'Concentration (IU/mL)': 'result_concentration',
        'Lower Limit of Detection (IU/mL)': 'llod',
        'Upper Limit of Detection (IU/mL)': 'ulod',
        'Testing Antigen': 'antigen',
        'Plate Reader': 'instrument',
        'Serial No.': 'instrument_serialno',
        'Software': 'lab_software_version'
    }
    df = df.rename(columns=rename)
    df = df.drop(columns=['Assay Lab','DER', 'instrument', 'lab_software_version'])

    ldms = pd.read_csv(constants.LDMS_PATH_HVTN,
                       dtype=constants.LDMS_DTYPE_MAP)
    ldms = ldms.loc[ldms.guspec.isin(df.guspec)]

    metadata_dict = {
            'network': 'HVTN',
            'upload_lab_id': 'P5',
            'assay_lab_name': 'Fouda Lab (Cornell)',
            'assay_type': 'ELISA',
            'assay_subtype': 'Hep B surface antigen',
            'specrole': 'Sample',
            'result_units': 'IU/mL',
            'lod_units': 'IU/mL',
            'instrument': 'BioTek Synergy LX Multimode Reader',
            'lab_software_version': 'BioTek Gen5 3.12',
    }

    outputs = sdmc.standard_processing(
            input_data=df,
            input_data_path=input_data_path,
            guspec_col="guspec",
            network="HVTN",
            metadata_dict=metadata_dict,
            ldms=ldms
        )

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
        'instrument_serialno',
        'lab_software_version',
        'lab_internal_sample_id',
        'assay_date',
        'plate_no',
        'antigen',
        'dilution',
        'result_concentration',
        'result_units',
        'llod',
        'ulod',
        'lod_units',
        'sdmc_processing_datetime',
        'sdmc_data_receipt_datetime',
        'input_file_name'
    ]

    mismatch = set(outputs.columns).symmetric_difference(reorder)
    if len(mismatch) > 0:
        raise Exception(f"trying to drop or add columns on reorder: {mismatch}")
    outputs = outputs[reorder]

    ## save to .txt
    savedir = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN135/assays/ELISA/misc_files/data_processing/'
    today = datetime.date.today().isoformat()
    fname = f"HVTN135_ELISA_HepB_Fouda_Processed_{today}.txt"

    outputs.to_csv(savedir + fname, sep="\t", index=False)

    # pivot summary
    pt1 = pd.pivot_table(outputs, index=['ptid', 'visitno'], columns='antigen', aggfunc='count')[['result_concentration']].droplevel(level=0, axis=1)
    pt1.to_excel(savedir + "HVTN135_HepB_ELISA_ptid_visitno_summary.xlsx")

    # check for completeness
    compare = pd.read_csv('/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN135/assays/bama_pvma/misc_files/data_processing/HVTN135_PVMA_FOUDA_processed_2024-08-27.txt', sep="\t")

    # exact match on guspecs
    set(compare.guspec).symmetric_difference(outputs.guspec)

    # the lab-internal id isn't consistent
    tmp = outputs[['guspec','lab_internal_sample_id']].merge(compare[['guspec','lab_internal_sample_id']], on=['guspec'])
    tmp.lab_internal_sample_id_y.nunique()
    tmp.loc[tmp.lab_internal_sample_id_x!=tmp.lab_internal_sample_id_y]

if __name__=="__main__":
    main()
