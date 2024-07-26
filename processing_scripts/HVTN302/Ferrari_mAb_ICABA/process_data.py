## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 07/26/2024
# Purpose: Process ICABA data from Ferrari lab
## ---------------------------------------------------------------------------##
import pandas as pd
import numpy as np
import datetime
import os

import sdmc_tools.process as sdmc
import sdmc_tools.constants as constants

def main():
    input_data_path = "/trials/vaccine/p302/s001/qdata/LabData/ICABA_mAb_pass-through/Ferrari_HVTN 302_mAbs_ICABA_Analysis_19JUL2024.csv"
    input_metadata_path = "/trials/vaccine/p302/s001/qdata/LabData/ICABA_mAb_pass-through/Ferrari_HVTN 302__mAbs_ICABA_Metadata_19JUL2024.csv"

    input_data = pd.read_csv(input_data_path)
    lab_metadata = pd.read_csv(input_metadata_path)

    lab_metadata["Antigen"] = lab_metadata.IMC_prep.str.replace(" ","_").str.split("_", expand=True)[0]

    data = data.rename(columns={
        'sample_id': 'mab_name',
        'Mock%IgG+': 'result_Mock%IgG+',
        '%p24+IgG+': 'result_%p24+IgG+',
        'MockSubtracted %p24+IgG+': 'result_MockSubtracted %p24+IgG+',
        'operator': 'lab_operator',
        'xml_file_name': 'lab_xml_file_name'
    })

    data.columns = [i.lower().replace(" ","_") for i in data.columns]

    # check the background subtraction done correctly
    # data['check'] = data["result_%p24+igg+"] -data['result_mock%igg+']
    # data.check[data.check < 0] = 0
    # data['error'] = np.abs(data.check-data['result_mocksubtracted_%p24+igg+'])
    # print(data.loc[data.error > 1e-10])

    md = {
        'network': ['HVTN'],
        'upload_lab_id': ['GF'],
        'assay_lab_name': ['Ferrari Lab'],
        'instrument': ['Fortessa Flow Cytometer'],
        'assay_type': ['ICABA'],
        'specrole': ['Monoclonal Antibody'],
        'result_units': ['Percent'],
    }
    data = data.merge(pd.DataFrame(md), how='cross')

    reorder = [
        'network',
        'protocol',
        'specrole',
        'mab_name',
        'upload_lab_id',
        'assay_lab_name',
        'assay_type',
        'instrument',
        'antigen',
        'concentration',
        'isotype',
        'result_mock%igg+',
        'result_%p24+igg+',
        'result_mocksubtracted_%p24+igg+',
        'result_units',
        'virus_stock',
        'virus_id',
        'imc_prep',
        'assay_date',
        'lab_operator',
        'fcs_file_name',
        'lab_xml_file_name',
        'analysis_file_name',
        'input_file_name',
        'input_metadata_file_name',
        'sdmc_processing_datetime',
        'sdmc_data_receipt_datetime',
    ]

    mismatch = set(data.columns).symmetric_difference(reorder)
    if len(mismatch) > 0:
        raise Exception("trying to add or remove columns!")

    outputs = data[reorder]
    outputs = outputs.sort_values(by=['mab_name', 'antigen', 'concentration'])

    savedir = "/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN302/assays/mAb_ICABA/misc_files/data_processing/"
    today = datetime.date.today().isoformat()
    save_fname = f"DRAFT_HVTN302_Ferrari_mAb_ICABA_processed_{today}.txt"

    outputs.to_csv(savedir + save_fname, index=False, sep="\t")

def pivot():
    summary = pd.pivot_table(data=outputs,
                   values='result_mocksubtracted_%p24+igg+',
                   index='mab_name',
                   columns=['antigen','concentration'],
                   aggfunc='count')
    summary.to_excel(savedir + "HVTN302_Ferrari_mab_icaba_pivot_summary.xlsx")

if __name__=="__main__":
    main()
    # pivot()
