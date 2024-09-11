## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 09/11/2024
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
    # read in ldms -----------------------------------------------------------##
    ldms = pd.read_csv(
        '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN115/specimens/ldms_feed/hvtn.ldms115.20240911.csv',
        usecols=constants.STANDARD_COLS,
        dtype=constants.LDMS_DTYPE_MAP
    )
    ldms = ldms.loc[ldms.lstudy==115.]

    # metadata ---------------------------------------------------------------##
    md = {
        'network': 'HVTN',
        'upload_lab_id': 'GF',
        'assay_lab_name': 'Ferrari Lab',
        'instrument': 'Fortessa Flow Cytometer',
        'assay_type': 'ICABA',
        'specrole': 'Sample',
        'result_units': 'Percent',
    }

    # process Part A (August 2024 data) --------------------------------------##
    # input data path
    input_data_path = '/trials/vaccine/p115/s001/qdata/Ferrari_ICABA_pass-through/Ferrari_HVTN 115_Part A_ICABA_Analysis_26AUG2024.csv'
    input_metadata_path = '/trials/vaccine/p115/s001/qdata/Ferrari_ICABA_pass-through/Ferrari_HVTN 115_Part A_ICABA_Metadata_26AUG2024.csv'

    input_data = pd.read_csv(input_data_path)
    input_metadata = pd.read_csv(input_metadata_path)

    input_data = input_data.dropna(how='all')
    input_metadata = input_metadata.dropna(how='all')

    partA = input_data.merge(input_metadata, on=['sample_id', 'ptid', 'visit', 'dilution'], how='outer')

    partA = partA.rename(columns={
        '%p24+IgG+': 'Infected %IgG+ Background Subtracted',
        'Mock%IgG+': 'Mock %IgG+ Background Subtracted',
        'MockSubtracted %p24+IgG+': 'Mock Subtracted %IgG+',
    })

    partA.columns = [i.lower().replace(" ","_") for i in partA.columns]

    partA = partA.rename(columns={
        'sample_id': 'guspec',
        'mock_%igg+_background_subtracted': 'result_mock_%igg+_background_subtracted',
        'infected_%igg+_background_subtracted': 'result_infected_%igg+_background_subtracted',
        'mock_subtracted_%igg+': 'result_mock_subtracted_%igg+',
        'xml_file_name': 'lab_xml_file_name',
    })

    # checked that these match ldms
    partA = partA.drop(columns=['ptid', 'visit'])

    partA['antigen'] = partA.imc_prep.str.split(" ", expand=True)[0]

    partA_outputs = sdmc.standard_processing(
        input_data=partA,
        input_data_path=input_data_path,
        guspec_col="guspec",
        network="HVTN",
        metadata_dict=md,
        ldms=ldms,
        additional_input_paths={'input_metadata_file_name': input_metadata_path}
    )

    # process Part B (Sep 2022 data) -----------------------------------------##
    # input data path
    partB_data_path = '/trials/vaccine/p115/s001/qdata/Ferrari_ICABA_pass-through/Ferrari_HVTN 115_ICABA_Analysis_20220921.csv'
    partB_metadata_path = '/trials/vaccine/p115/s001/qdata/Ferrari_ICABA_pass-through/Ferrari_HVTN 115_ICABA_Metadata_20220921.csv'

    partB_data = pd.read_csv(partB_data_path)
    partB_metadata = pd.read_csv(partB_metadata_path)

    partB_data = partB_data.dropna(how='all')

    partB = partB_data.merge(partB_metadata, on=['sample_id', 'ptid', 'visit', 'dilution'], how='outer')

    #checked that these match ldms
    partB = partB.drop(columns=['visit','ptid'])
    partB.columns = [i.lower().replace(" ","_") for i in partB.columns]

    partB = partB.rename(columns={
        'sample_id': 'guspec',
        'mock_%igg+_background_subtracted': 'result_mock_%igg+_background_subtracted',
        'infected_%igg+_background_subtracted': 'result_infected_%igg+_background_subtracted',
        'mock_subtracted_%igg+': 'result_mock_subtracted_%igg+',
        'cd4+_%igg+_background_subtracted': 'result_cd4+_%igg+_background_subtracted',
        'cd4-_%igg+_background_subtracted': 'result_cd4-_%igg+_background_subtracted',
        'xml_file_name': 'lab_xml_file_name',
    })

    partB_outputs = sdmc.standard_processing(
        input_data=partB,
        input_data_path=partB_data_path,
        guspec_col="guspec",
        network="HVTN",
        metadata_dict=md,
        ldms=ldms,
        additional_input_paths={'input_metadata_file_name': partB_metadata_path}
    )

    # check mock subtraction
    mock_subtracted = outputs['result_infected_%igg+_background_subtracted'] - outputs['result_mock_%igg+_background_subtracted']
    mock_subtracted[mock_subtracted < 0] = 0
    discrep = np.abs(outputs['result_mock_subtracted_%igg+'] - mock_subtracted)
    if discrep.max() > 1e-10:
        raise Exception("Mock subtraction doesn't look correct")

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
        'result_infected_%igg+_background_subtracted',
        'result_mock_%igg+_background_subtracted',
        'result_mock_subtracted_%igg+',
        'result_cd4+_%igg+_background_subtracted',
        'result_cd4-_%igg+_background_subtracted',
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

    savedir = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN115/assays/icaba/misc_files/data_processing/'
    today = datetime.date.today().isoformat()
    outputs.to_csv(savedir + f"DRAFT_HVTN115_ICABA_FERRARI_processed_{today}.txt", sep="\t", index=False)

    pivot_summary = pd.pivot_table(data=outputs, index='ptid', columns='visitno', aggfunc='count', fill_value=0)['result_mock_subtracted_%igg+']
    pivot_summary.to_excel(savedir + "HVTN115_ICABA_ptid_visitno_summary.xlsx")

    partwise_pivot_summary = pd.pivot_table(data=outputs, index='ptid', columns=['input_file_name','visitno'], aggfunc='count', fill_value=0)['result_mock_subtracted_%igg+']
    partwise_pivot_summary.to_excel(savedir + "HVTN115_ICABA_part_ptid_visitno_summary.xlsx")


if __name__=="__main__":
    main()
