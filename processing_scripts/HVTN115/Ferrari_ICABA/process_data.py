## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 08/26/2024
# Purpose:  Process ICABA data from Ferrari lab for HVTN115
## ---------------------------------------------------------------------------##
import pandas as pd
import numpy as np
import datetime
import yaml
import os

import sdmc_tools.refactor_process as sdmc
import sdmc_tools.constants as constants
## ---------------------------------------------------------------------------##

def main():
    # input data path
    input_data_path = '/trials/vaccine/p115/s001/qdata/Ferrari_ICABA_pass-through/Ferrari_HVTN 115_Part A_ICABA_Analysis_26AUG2024.csv'
    input_metadata_path = '/trials/vaccine/p115/s001/qdata/Ferrari_ICABA_pass-through/Ferrari_HVTN 115_Part A_ICABA_Metadata_26AUG2024.csv'

    input_data = pd.read_csv(input_data_path)
    input_metadata = pd.read_csv(input_metadata_path)

    input_data = input_data.dropna(how='all')
    input_metadata = input_metadata.dropna(how='all')

    data = input_data.merge(input_metadata, on=['sample_id', 'ptid', 'visit', 'dilution'], how='outer')
    data.columns = [i.lower().replace(" ","_") for i in data.columns]

    # check mock subtraction
    mock_subtracted = data['%p24+igg+'] - data['mock%igg+']
    mock_subtracted[mock_subtracted < 0] = 0
    discrep = sum(data['mocksubtracted_%p24+igg+'] - mock_subtracted)
    if discrep > 0:
        raise Exception("Mock subtraction doesn't look correct")

    data = data.rename(columns={
        'sample_id': 'guspec',
        'mock%igg+': 'result_mock%igg+',
        '%p24+igg+': 'result_%p24+igg+',
        'mocksubtracted_%p24+igg+': 'result_mocksubtracted_%p24+igg+',
        'xml_file_name': 'lab_xml_file_name',
    })

    # checked that these match ldms
    data = data.drop(columns=['ptid', 'visit'])

    data['antigen'] = data.imc_prep.str.split(" ", expand=True)[0]

    ldms = pd.read_csv(
        constants.LDMS_PATH_HVTN,
        usecols=constants.STANDARD_COLS,
        dtype=constants.LDMS_DTYPE_MAP
    )
    ldms = ldms.loc[ldms.lstudy==115.]
    ldms = ldms.loc[ldms.guspec.isin(data.guspec)]

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
        guspec_col="guspec",
        network="HVTN",
        metadata_dict=md,
        ldms=ldms,
        additional_input_paths={'input_metadata_file_name': input_metadata_path}
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
        'instrument',
        'antigen',
        'isotype',
        'dilution',
        'result_%p24+igg+',
        'result_mock%igg+',
        'result_mocksubtracted_%p24+igg+',
        'result_units',
        'assay_date',
        'operator',
        'virus_stock',
        'imc_prep',
        'fcs_file_name',
        'analysis_file_name',
        'lab_xml_file_name',
        'sdmc_processing_datetime',
        'sdmc_data_receipt_datetime',
        'input_file_name',
        'input_metadata_file_name',
    ]
    mismatch = set(outputs.columns).symmetric_difference(reorder)

    if len(mismatch) > 0:
        raise Exception("Trying to add or remove columns")
    outputs = outputs[reorder]

    ## save to txt
    savedir = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN115/assays/icaba/misc_files/data_processing/'
    today = datetime.date.today().isoformat()
    outputs.to_csv(savedir + f"DRAFT_HVTN115_ICABA_FERRARI_processed_{today}.txt", sep="\t", index=False)

    ## save pivot summary
    pivot_summary = pd.pivot_table(data=outputs, index='ptid', columns='visitno', aggfunc='count', fill_value=0)['result_mocksubtracted_%p24+igg+']
    pivot_summary.to_excel(savedir + "HVTN115_ICABA_ptid_visitno_summary.xlsx")


if __name__=="__main__":
    main()
