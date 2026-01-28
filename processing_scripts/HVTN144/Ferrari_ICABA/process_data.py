## ---------------------------------------------------------------------------##
# Author: Sara Thiebaud
# Date: 22 December 2025
# Purpose:  Process ICABA data from Ferrari lab for HVTN144
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
    # read in ldms -----------------------------------------------------------##
    ldms = pd.read_csv(
        constants.LDMS_PATH_HVTN,
        usecols=constants.STANDARD_COLS,
        dtype=constants.LDMS_DTYPE_MAP
    )
    ldms = ldms.loc[ldms.lstudy==144.]

    # input data paths -------------------------------------------------------##
    metadata_path = '/trials/vaccine/p144/s001/qdata/LabData/ICABA_pass-through/20251219-02/Ferrari_HVTN 144_ICABA_Metadata_19DEC2025.csv'
    metadata = pd.read_csv(metadata_path)

    input_data_path = '/trials/vaccine/p144/s001/qdata/LabData/ICABA_pass-through/20251219-02/Ferrari_HVTN 144_ICABA_Analysis_19DEC2025.csv'
    data = pd.read_csv(input_data_path)

    data = data.rename(columns={"sample_id":"guspec",
                            "visit":"visitno",
                            'Mock%IgG+': 'result_Mock%IgG+',
                            '%p24+IgG+': 'result_%p24+IgG+',
                            'MockSubtracted %p24+IgG+': 'result_MockSubtracted %p24+IgG+'
                           })

    metadata = metadata.rename(columns={'sample_id':'guspec',
                                    'visit':'visitno'})

    data = data.merge(metadata,
                      on=["guspec", "ptid", "visitno", "dilution"],
                      how='outer'
    )

    # going to make sure these are the same
    data = data.rename(columns={
        'ptid':'ptid_compare',
        'visitno':'visitno_compare'
    })

    # reformat data + std processing -----------------------------------------##
    data = data.drop(columns=['analysis_file_name', 'result_MockSubtracted %p24+IgG+'])
    data = data.rename(columns={'xml_file_name':'lab_xml_file_name'})

    md = {
        'network': 'HVTN',
        'upload_lab_id': 'GF',
        'assay_lab_name': 'Ferrari Lab',
        'instrument': 'Fortessa Flow Cytometer',
        'assay_type': 'ICABA',
        'specrole': 'Sample',
        'result_units': 'Percent',
    }


    ldms = ldms.loc[ldms.guspec.isin(data.guspec)]
    outputs = sdmc.standard_processing(
        input_data=data,
        input_data_path=input_data_path,
        additional_input_paths={'input_metadata_file_name':metadata_path},
        guspec_col="guspec",
        network="HVTN",
        metadata_dict=md,
        ldms=ldms
    )

    # ensure these are the same, then drop
    outputs.visitno = outputs.visitno.astype(float)
    outputs.visitno_compare = outputs.visitno_compare.astype(float)
    if (outputs.visitno!=outputs.visitno_compare).sum() > 0:
        raise Warning("lab submitted visitnos don't match ldms")

    outputs.ptid = outputs.ptid.astype(int)
    outputs.ptid_compare = outputs.ptid_compare.astype(int)
    if (outputs.ptid!=outputs.ptid_compare).sum() > 0:
        raise Warning("lab submitted ptids don't match ldms")

    outputs = outputs.drop(columns=['ptid_compare','visitno_compare'])

    # final formatting -------------------------------------------------------##
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
        'instrument',
        'dilution',
        'result_%p24+igg+',
        'result_mock%igg+',
        #'result_mocksubtracted_%p24+igg+',
        'result_units',
        'virus_id',
        'virus_stock',
        'imc_prep',
        'assay_date',
        'fcs_file_name',
        'isotype',
        'lab_xml_file_name',
        'operator',
        'sdmc_processing_datetime',
        'sdmc_data_receipt_datetime',
        'input_file_name',
        'input_metadata_file_name'
    ]

    mismatch = set(outputs.columns).symmetric_difference(reorder)

    if len(mismatch) > 0:
        raise Exception(f"trying to drop or add columns on reorder: {mismatch}")
    outputs = outputs[reorder]

    check = outputs['result_%p24+igg+'] - outputs['result_mock%igg+']
    check[check < 0] = 0

    # if np.abs(check - outputs['result_mocksubtracted_%p24+igg+']).max() > 1e-10:
    #     raise Exception("Mock subtraction looks wrong")

    # save results -----------------------------------------------------------##
    savedir = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN144/assays/icaba/misc_files/data_processing/'
    today = datetime.date.today().isoformat()
    fname = f"HVTN144_Ferrari_ICABA_Processed_{today}.txt"

    outputs.to_csv(savedir + fname, sep="\t", index=False)

    # pivot summary ----------------------------------------------------------##
    summary = pd.pivot_table(outputs, index='ptid', columns='visitno', aggfunc='count', fill_value=0)['result_mock%igg+']
    summary.to_excel(savedir + "HVTN144_ICABA_summary.xlsx")

if __name__=="__main__":
    main()
