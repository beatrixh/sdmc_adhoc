## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 06/10/2024
# Purpose:  - Process SECABA data from Ferrari lab
## ---------------------------------------------------------------------------##
import pandas as pd
import datetime, os
import yaml
import sdmc_tools.process as sdmc
import sdmc_tools.constants as constants

def main():
    ## read in data ----------------------------------------------------------##
    secaba_analysis_path = "/trials/covpn/p3008/s001/qdata/LabData/SECABA_pass-through/Ferrari_CoVPN3008_SECABA_Analysis_2024-01-02.csv"
    secaba_metadata_path = "/trials/covpn/p3008/s001/qdata/LabData/SECABA_pass-through/Ferrari_CoVPN3008_SECABA_Metadata_2024-01-25.csv"

    data = pd.read_csv(secaba_analysis_path)
    metadata = pd.read_csv(secaba_metadata_path)

    data.columns = [i.lower().replace(" ", "_").replace("%", "pct_") for i in data.columns]

    rename = {'sample_id':'guspec', 'ptid':'submitted_ptid', 'visit': 'submitted_visit'}
    data = data.rename(columns=rename)
    metadata = metadata.rename(columns=rename)

    # drop this for now until have clarification from lab
    data = data.loc[data.guspec!="0392-01JWBJ00-002"]
    metadata = metadata.loc[metadata.guspec!="0392-01JWBJ00-002"]

    # merge metadata onto data
    data = data.merge(metadata, on=['guspec', 'submitted_ptid', 'submitted_visit', 'dilution'], how='outer')

    ## hand-entered metadata
    metadata_dict = {
        'network': 'CoVPN',
        'upload_lab_id': 'GF',
        'assay_lab_name': 'Ferrari Lab',
        'instrument': 'Fortessa Flow Cytometer',
        'assay_type': 'SECABA',
        'specrole': 'Sample',
        'result_units': 'Percent',
    }

    ldms = pd.read_csv(
        constants.LDMS_PATH_COVPN,
        usecols=constants.STANDARD_COLS,
        dtype=constants.LDMS_DTYPE_MAP
    )
    ldms = ldms.loc[ldms.lstudy==3008.]
    ldms = ldms.loc[ldms.guspec.isin(data.guspec)]

    ldms.txtpid = ldms.txtpid.astype(str).astype(int)
    ldms = ldms.drop_duplicates()

    outputs = sdmc.standard_processing(
        input_data=data,
        input_data_path=secaba_analysis_path,
        guspec_col='guspec',
        network='covpn',
        metadata_dict=metadata_dict,
        ldms=ldms
    )
    outputs = outputs.drop(columns=["submitted_visit", "submitted_ptid"])

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
        'mock_pct_igg+',
        'transfected_pct_igg+',
        'background_subtracted_pct_igg+',
        'result_units',
        'isotype',
        'assay_date',
        'operator',
        'analysis_file_name',
        'fcs_file_name',
        'xml_file_name',
        'sdmc_processing_datetime',
        'sdmc_data_receipt_datetime',
        'input_file_name'
    ]

    # set(reorder).symmetric_difference(outputs.columns)
    outputs = outputs[reorder]
    outputs['input_file_name'] = "Ferrari_CoVPN3008_SECABA_Analysis_2024-01-02.csv, SECABA_pass-through/Ferrari_CoVPN3008_SECABA_Metadata_2024-01-25.csv"

    savedir = "/networks/vtn/lab/SDMC_labscience/studies/CoVPN/CoVPN3008/assays/SECABA/misc_files/data_processing/"
    today = datetime.datetime.today().strftime('%Y-%m-%d')
    fname = f"CoVPN3008_Ferrari_SECABA_processed_{today}.txt"
    outputs.to_csv(savedir + fname, index=False, sep="\t")


if __name__ == '__main__':
    main()
