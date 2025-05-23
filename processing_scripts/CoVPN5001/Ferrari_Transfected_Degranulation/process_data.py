## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 11-14-2024
# Purpose:  - Process CoVPN5001 Transfected Degranulation data from Ferrari lab
## ---------------------------------------------------------------------------##
import pandas as pd
import datetime, os
import numpy as np
import yaml
import sdmc_tools.process as sdmc
import sdmc_tools.constants as constants
## ---------------------------------------------------------------------------##

def main():
    transfected_degran_analysis_path = "/trials/covpn/p5001/s001/qdata/LabData/Ferrari_pass-through/Ferrari_CoVPN5001_Transfected Cell Degranulation_Analysis_2022-11-15.csv"
    transfected_degran_metadata_path = "/trials/covpn/p5001/s001/qdata/LabData/Ferrari_pass-through/Ferrari_CoVPN5001_Transfected Cell Degranulation_Metadata_2023-01-17.csv"

    data = pd.read_csv(transfected_degran_analysis_path)
    metadata = pd.read_csv(transfected_degran_metadata_path)

    data = data.merge(metadata, on=['sample_id','ptid','visit','dilution'], how='outer')
    data.columns = [i.lower().replace(" ", "_").replace("%", "pct_") for i in data.columns]

    data = data.rename(columns={
        'sample_id': 'guspec',
        'ptid': 'ptid_compare',
        'visit': 'visit_compare',
        'xml_file_name': 'lab_xml_file_name',
        'infected_pct_cd107a+_nk_cells': 'result_infected_pct_cd107a+_nk_cells',
        'mock_pct_cd107a+_nk_cells': 'result_mock_pct_cd107a+_nk_cells',
        'background_subtracted_pct_cd107a+_nk_cells': 'result_background_subtracted_pct_cd107a+_nk_cells',
    })

    ldms = pd.read_csv(
        constants.LDMS_PATH_HVTN, ##yes this is correct, idk why this is here
        usecols=constants.STANDARD_COLS,
        dtype=constants.LDMS_DTYPE_MAP
    )
    ldms = ldms.loc[ldms.guspec.isin(data.guspec)]

    ldms.txtpid = ldms.txtpid.astype(str).astype(int)
    ldms = ldms.drop_duplicates()

    ## hand-entered metadata
    metadata_dict = {
        'network': 'CoVPN',
        'upload_lab_id': 'GF',
        'assay_lab_name': 'Ferrari Lab',
        'instrument': 'Fortessa Flow Cytometer',
        'assay_type': 'Antibody-Dependent Cellular Cytotoxicity (ADCC)',
        'assay_subtype': 'Transfected Cell CD107a Degranulation',
        'assay_details': 'Transfected',
        'specrole': 'Sample',
        'result_units': 'Percent',
    }

    outputs = sdmc.standard_processing(
        input_data=data,
        input_data_path=transfected_degran_analysis_path,
        guspec_col='guspec',
        network='covpn',
        metadata_dict=metadata_dict,
        ldms=ldms,
        additional_input_paths={'input_metadata_file_name': transfected_degran_metadata_path}
    )

    outputs.ptid = outputs.ptid.astype(str)
    outputs.ptid_compare = outputs.ptid_compare.astype(str)

    ptid_diff = sum(outputs.ptid != outputs.ptid_compare)
    if ptid_diff > 0:
        raise Warning("Lab-submitted ptid doesn't match ldms")

    outputs.visitno = outputs.visitno.astype(float)
    outputs.visit_compare = outputs.visit_compare.astype(float)

    visit_diff = sum(outputs.visitno != outputs.visit_compare)
    if visit_diff > 0:
        raise Warning("Lab-submitted visitno doesn't match ldms")

    outputs = outputs.drop(columns=['visit_compare','ptid_compare','infected/transfected'])

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
        'dilution',
        'result_infected_pct_cd107a+_nk_cells',
        'result_mock_pct_cd107a+_nk_cells',
        'result_background_subtracted_pct_cd107a+_nk_cells',
        'result_units',
        'assay_date',
        'operator',
        'analysis_file_name',
        'fcs_file_name',
        'lab_xml_file_name',
        'sdmc_processing_datetime',
        'sdmc_data_receipt_datetime',
        'input_file_name',
        'input_metadata_file_name',
    ]

    mismatch = set(reorder).symmetric_difference(outputs.columns)
    if len(mismatch) > 0:
        raise Exception("Trying to add or remove columns")
    outputs = outputs[reorder]

    savedir = "/networks/vtn/lab/SDMC_labscience/studies/CoVPN/CoVPN5001/assays/ADCC/misc_files/data_processing/transfected_degranulation/"
    today = datetime.datetime.today().strftime('%Y-%m-%d')
    fname = f"CoVPN5001_Ferrari_Transfected_Degranulation_processed_{today}.txt"
    outputs.to_csv(savedir + fname, index=False, sep="\t")

    summary = pd.pivot_table(outputs,
                             index=['ptid'],
                             columns=['dilution','visitno'],
                             aggfunc='count',
                             fill_value=0)['result_background_subtracted_pct_cd107a+_nk_cells']
    summary.to_excel(savedir + "CoVPN5001_Transfected_Degran_summary.xlsx")

if __name__=="__main__":
    main()
