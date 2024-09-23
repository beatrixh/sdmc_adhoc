## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 09/23/2024
# Purpose:  - Process ICABA data from Ferrari Lab
## ---------------------------------------------------------------------------##
import pandas as pd
import numpy as np
import datetime
import yaml
import os

import sdmc_tools.refactor_process as sdmc
import sdmc_tools.constants as constants

def main():
    # input data paths -------------------------------------------------------##
    input_data_path = '/trials/vaccine/p305/s001/qdata/LabData/ICABA_pass-through/Ferrari_HVTN 305_ICABA_Analysis_08AUG2024.csv'
    metadata_path = '/trials/vaccine/p305/s001/qdata/LabData/ICABA_pass-through/Ferrari_HVTN 305_ICABA_Metadata_08AUG2024.csv'

    # read in ldms -----------------------------------------------------------##
    ldms = pd.read_csv(
        constants.LDMS_PATH_HVTN,
        usecols=constants.STANDARD_COLS,
        dtype=constants.LDMS_DTYPE_MAP
    )
    ldms = ldms.loc[ldms.lstudy==305.]

    # reformat data + std processing -----------------------------------------##
    data = pd.read_csv(input_data_path)
    metadata = pd.read_csv(metadata_path)

    data = data.rename(columns={"sample_id":"guspec",
                            "visit":"visitno",
                            'Mock%IgG+': 'result_mock_%igg+_background_subtracted',
                            '%p24+IgG+': 'result_infected_%igg+_background_subtracted',
                            'MockSubtracted %p24+IgG+': 'result_mock_subtracted_%igg+'
                           })
    metadata = metadata.rename(columns={'sample_id':'guspec',
                                    'visit':'visitno'})
    data = data.merge(metadata,
                      on=["guspec", "ptid", "visitno", "dilution"],
                      how='outer'
    )

    # checked that these match ldms
    data = data.drop(columns=["visitno", "ptid", "analysis_file_name"])
    data = data.rename(columns={'xml_file_name':'lab_xml_file_name'})

    data['antigen'] = data.IMC_prep.str.split("_", expand=True)[0]

    ldms = ldms.loc[ldms.guspec.isin(data.guspec)].drop_duplicates()

    md = {
        'network': 'HVTN',
        'upload_lab_id': 'GF',
        'assay_lab_name': 'Ferrari Lab',
        'instrument': 'Fortessa Flow Cytometer',
        'assay_type': 'ICABA',
        'specrole': 'Sample',
        'result_units': 'Percent',
    }

    outputs = sdmc.standard_processing(
        input_data=data,
        input_data_path=input_data_path,
        additional_input_paths={'input_metadata_file_name':metadata_path},
        guspec_col="guspec",
        network="HVTN",
        metadata_dict=md,
        ldms=ldms
    )

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
        'antigen',
        'dilution',
        'result_infected_%igg+_background_subtracted',
        'result_mock_%igg+_background_subtracted',
        'result_mock_subtracted_%igg+',
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

    check = outputs['result_infected_%igg+_background_subtracted'] - outputs['result_mock_%igg+_background_subtracted']
    check[check < 0] = 0

    if np.abs(check - outputs['result_mock_subtracted_%igg+']).max() > 1e-10:
        raise Exception("Mock subtraction looks wrong")

    # save results -----------------------------------------------------------##
    savedir = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN305/assays/ICABA/misc_files/data_processing/'
    today = datetime.date.today().isoformat()
    fname = f"HVTN305_Ferrari_ICABA_Processed_{today}.txt"

    outputs.to_csv(savedir + fname, sep="\t", index=False)

    # pivot summary ----------------------------------------------------------##
    summary = pd.pivot_table(outputs, index=['ptid'], columns='visitno', aggfunc='count', fill_value=0)['result_mocksubtracted_%p24+igg+']
    summary.to_excel(savedir + "HVTN305_ICABA_summary.xlsx")

    # manifest checks --------------------------------------------------------##

    m1 = pd.read_csv("/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN305/assays/ICABA/misc_files/hvtn305_ferrari_manifest1.txt", sep="\t")
    m2 = pd.read_csv("/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN305/assays/ICABA/misc_files/hvtn305_ferrari_manifest2.txt", sep="\t")
    ms = pd.concat([m1, m2])

    manids = set(m1.GLOBAL_ID).union(m2.GLOBAL_ID)

    mismatch = manids.symmetric_difference(outputs.guspec)
    if len(mismatch) > 0:
        print(f"Manifest/data don't match up")
        print(mismatch)

        print(ms.loc[ms.GLOBAL_ID.isin(mismatch),['GLOBAL_ID','PROTOCOL','PID','VID']])

if __name__=="__main__":
    main()
