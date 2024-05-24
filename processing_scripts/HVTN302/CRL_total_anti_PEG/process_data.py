## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 02/29/2024
# Purpose:  - Process CRL total anti-PEG data processing
## ---------------------------------------------------------------------------##
import pandas as pd
import datetime, os
import yaml
import sdmc_tools.process as sdmc
import sdmc_tools.constants as constants

## Read in data --------------------------------------------------------------##
def main():
    ldms = pd.read_csv(constants.LDMS_PATH_HVTN)
    ldms = ldms.loc[ldms.lstudy==302.]

    yaml_path = os.path.dirname(__file__) + "/paths.yaml"
    with open(yaml_path, 'r') as file:
        yaml_dict = yaml.safe_load(file)

    input_data_path = yaml_dict['input_data_path']
    data = pd.read_excel(input_data_path)

    # add guspec column
    data['guspec'] = data['Sample \nUser text 1'].str[3:]

    # standardize column format
    data.columns = [i.lower().replace("\n", "").replace(" ", "_") for i in data.columns]

    # map results column
    titer_map = {
        100: '100',
        400: '400',
        1600: '1600',
        'Positive Titer (>6400.000)': '>6400',
        3200: '3200',
        'Negative Screen': "N/A (Negative Screen)",
        'Negative Titer (<50.000)': '<50',
        800: '800',
        200: '200',
        'Negative Immunodepletion': 'N/A (Negative Immunodepletion)',
        50: '50'
    }

    interpretation_map = {
        100: 'Positive',
        400: 'Positive',
        1600: 'Positive',
        'Positive Titer (>6400.000)': 'Positive',
        3200: 'Positive',
        'Negative Screen': 'Negative',
        'Negative Titer (<50.000)': 'Negative',
        800: 'Positive',
        200: 'Positive',
        'Negative Immunodepletion': 'Negative',
        50: 'Positive'
    }

    units_map = {
        100: 'Titer',
        400: 'Titer',
        1600: 'Titer',
        'Positive Titer (>6400.000)': 'Titer',
        3200: 'Titer',
        'Negative Screen': 'N/A (Negative Screen)',
        'Negative Titer (<50.000)': 'Titer',
        800: 'Titer',
        200: 'Titer',
        'Negative Immunodepletion': 'N/A (Negative Immunodepletion)',
        50: 'Titer'
    }

    data["result"] = data.ir_result_text.map(titer_map)
    data["result_interpretation"] = data.ir_result_text.map(interpretation_map)
    data["result_units"] = data.ir_result_text.map(units_map)

    ## standard processing
    metadata_dict = {
        'network': 'HVTN',
        'upload_lab_id': 'N4',
        'assay_lab_name': 'Charles River Laboratories (CRL)',
        'instrument': 'SpectraMax',
        'assay_type': 'Total Anti-PEG',
        'specrole': 'Sample',
        'assay_precision': 'Qualitative',
    }

    outputs = sdmc.standard_processing(
        input_data = data,
        input_data_path=input_data_path,
        guspec_col='guspec',
        network='hvtn',
        metadata_dict=metadata_dict,
        ldms=ldms
    )

    ## checks ##

    if sum(outputs.visit.astype(int) != outputs.visitno.astype(float).astype(int)) > 0:
        raise Exception("visitnos don't match ldms")

    if sum(outputs.subject.astype(int)!=outputs.ptid.astype(int)) > 0:
        raise Exception("ptids don't match ldms")

    ## drop columns
    outputs = outputs.drop(columns=[
        "subject",
        "visit",
        'day_nominal',
        'nominal_time',
        'study_day',
        'study_id',
        'project_name',
        'barcode_id',
        'custom_id',
        'sample_name',
        'ir_result_text'
    ])

    ## rename columns
    rename = {
        'run_id': 'lab_internal_run_id',
        'sample_receipt_date': 'lab_internal_sample_receipt_date',
        'sample_user_text_1': 'sample_id_submitted',
    }

    outputs = outputs.rename(columns=rename)

    ## reorder
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
        'analyte',
        'result',
        'result_units',
        'result_interpretation',
        'assay_precision',
        'lab_internal_run_id',
        'lab_internal_sample_receipt_date',
        'sample_id_submitted',
        'sdmc_processing_datetime',
        'sdmc_data_receipt_datetime',
        'input_file_name'
    ]

    if len(list(dir(set(reorder).symmetric_difference(outputs.columns)))) > 0:
        raise Exception("Oops; trying to drop or add columns")

    outputs = outputs[reorder]

    ## save outputs
    savedir = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN302/assays/AE_assays/Anti-PEG_total/misc_files/data_processing/'
    today = datetime.datetime.today().strftime('%Y-%m-%d')
    fname = f"DRAFT_HVTN302_Anti-PEG_total_CRL_processed_{today}.txt"
    outputs.to_csv(savedir + fname, sep="\t", index=False)


def pivot():
    pd.pivot_table(data=outputs, values='result', index='ptid', columns='visitno', aggfunc='count', fill_value=0).to_csv(savedir + "ptid_visitno_summary.csv")
    outputs.guspec.value_counts().to_frame().to_csv(savedir + "guspec_summary.csv")