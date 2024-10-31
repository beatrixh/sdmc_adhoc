## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 10/25/2024
# Purpose: Process transfected degranulation data from Ferrari lab
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
    ldms = pd.read_csv(
        constants.LDMS_PATH_HVTN,
        usecols=constants.STANDARD_COLS,
        dtype=constants.LDMS_DTYPE_MAP
    )
    ldms = ldms.loc[ldms.lstudy==405.]

    input_data_path = '/trials/vaccine/p405/s001/qdata/Ferrari/Ferrari_HVTN405_Transfected_cell_degranulation_analysis_2021-09-15.csv'
    metadata_path = '/trials/vaccine/p405/s001/qdata/Ferrari/Ferrari_HVTN405_Transfected_cell_degranulation_metadata_2021-09-15.csv'

    metadata = pd.read_csv(metadata_path)
    data = pd.read_csv(input_data_path)

    data = data.merge(metadata, on=['sample_id','ptid','dilution'], how='outer')
    data.columns = [i.lower().replace(" ","_") for i in data.columns]

    data = data.rename(columns={
        'sample_id':'guspec',
        'mock_%cd107a+_nk_cells': 'result_mock_%cd107a+_nk_cells',
        'infected_%cd107a+_nk_cells': 'result_infected_%cd107a+_nk_cells',
        'background_subtracted_%cd107a+_nk_cells': 'result_background_subtracted_%cd107a+_nk_cells',
        'xml_file_name': 'lab_xml_file_name',
        'visit':'compare_visit',
        'ptid':'compare_ptid',
    })

    check = data['result_infected_%cd107a+_nk_cells'] - data['result_mock_%cd107a+_nk_cells']
    check[check < 0] = 0

    background_subtraction_check = np.abs(check - data['result_background_subtracted_%cd107a+_nk_cells']).max()
    if background_subtraction_check > 1e-10:
        raise Warning("Background subtraction looks incorrect")

    # metadata from sdmc background materials
    metadata_dict = {
        'network': 'HVTN',
        'upload_lab_id': 'GF',
        'assay_lab_name': 'Ferrari Lab',
        'instrument': 'Fortessa Flow Cytometer',
        'assay_type': 'Antibody-Dependent Cellular Cytotoxicity (ADCC)',
        'assay_subtype': 'Transfected Cell CD107a Degranulation',
        'assay_details': 'Transfected',
        'specrole': 'Sample',
        'result_units': 'Percent',
    }

    # standard processing
    outputs = sdmc.standard_processing(
        input_data=data,
        input_data_path=input_data_path,
        guspec_col='guspec',
        network='hvtn',
        metadata_dict=metadata_dict,
        ldms=ldms.loc[ldms.guspec.isin(data.guspec)]
    )

    outputs.compare_visit = outputs.compare_visit.astype(float)
    outputs.visitno = outputs.visitno.astype(float)
    diff = (outputs.visitno!=outputs.compare_visit).sum()

    if diff > 0:
        raise Warning("visitno doesn't match lab submission")

    outputs.compare_ptid = outputs.compare_ptid.astype(str)
    outputs.ptid = outputs.ptid.astype(str)
    diff = (outputs.compare_ptid!=outputs.ptid).sum()

    if diff > 0:
        raise Warning("ptid doesn't match lab submission")


    outputs = outputs.drop(columns=['compare_visit','compare_ptid'])
    outputs = outputs.drop(columns=['infected/transfected'])

    # reorder outputs
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
        'assay_details',
        'instrument',
        # 'antigen',
        'dilution',
        'result_infected_%cd107a+_nk_cells',
        'result_mock_%cd107a+_nk_cells',
        'result_background_subtracted_%cd107a+_nk_cells',
        'result_units',
        'assay_date',
        'operator',
        'analysis_file_name',
        'fcs_file_name',
        'lab_xml_file_name',
        'sdmc_processing_datetime',
        'sdmc_data_receipt_datetime',
        'input_file_name',
    ]

    mismatch = set(reorder).symmetric_difference(outputs.columns)
    if len(mismatch) > 0:
        raise Warning("Trying to add or remove columns")
    outputs = outputs[reorder]

    savedir = "/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN405_HPTN1901/assays/tranfected_cell_cd107a_degranulation/misc_files/data_processing/"
    today = datetime.datetime.today().strftime('%Y-%m-%d')
    fname = f"HVTN405_Ferrari_ADCC_Degranulation_processed_{today}.txt"
    outputs.to_csv(savedir + fname, sep="\t", index=False)

    summary = pd.pivot_table(outputs, index='ptid', columns='visitno', aggfunc='count', fill_value=0)['result_background_subtracted_%cd107a+_nk_cells']
    summary.to_excel(savedir + "HVTN405_data_summary.xlsx")

if __name__=="__main__":
    main()
