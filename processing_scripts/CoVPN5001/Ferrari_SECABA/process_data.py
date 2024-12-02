## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 10-30-2024
# Purpose:  - Process CoVPN5001 SECABA data from Ferrari lab
## ---------------------------------------------------------------------------##
import pandas as pd
import datetime, os
import yaml
import sdmc_tools.process as sdmc
import sdmc_tools.constants as constants
## ---------------------------------------------------------------------------##

def main():
    secaba_analysis_path = "/trials/covpn/p5001/s001/qdata/LabData/Ferrari_pass-through/Ferrari_CoVPN5001_Infected Cell Degranulation_Analysis_2022-03-23.csv"
    secaba_metadata_path = "/trials/covpn/p5001/s001/qdata/LabData/Ferrari_pass-through/Ferrari_CoVPN5001_Infected Cell Degranultaion_Metadata_2022-03-25.csv"

    data = pd.read_csv(secaba_analysis_path)
    metadata = pd.read_csv(secaba_metadata_path)

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
        constants.LDMS_PATH_HVTN,
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
        'assay_type': 'SECABA',
        'specrole': 'Sample',
        'result_units': 'Percent',
    }

    outputs = sdmc.standard_processing(
        input_data=data,
        input_data_path=secaba_analysis_path,
        guspec_col='guspec',
        network='covpn',
        metadata_dict=metadata_dict,
        ldms=ldms,
        additional_input_paths={'input_metadata_file_name': secaba_metadata_path}
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

    outputs = outputs.drop(columns=['visit_compare','ptid_compare'])

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
        'infected/transfected',
        'instrument',
        'dilution',
        'result_infected_pct_cd107a+_nk_cells',
        'result_mock_pct_cd107a+_nk_cells',
        'result_background_subtracted_pct_cd107a+_nk_cells',
        'result_units',
        # 'isotype',
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

    savedir = "/networks/vtn/lab/SDMC_labscience/studies/CoVPN/CoVPN5001/assays/ADCC/misc_files/data_processing/SECABA/"
    today = datetime.datetime.today().strftime('%Y-%m-%d')
    fname = f"CoVPN5001_Ferrari_SECABA_processed_{today}.txt"
    outputs.to_csv(savedir + fname, index=False, sep="\t")

    summary = pd.pivot_table(outputs, index=['ptid'], columns=['dilution','visitno'], aggfunc='count', fill_value=0)['result_background_subtracted_pct_cd107a+_nk_cells']
    summary.to_excel(savedir + "CoVPN5001_SECABA_data_summary.xlsx")

    ## check for completeness agasint manifest
    most_samples_path = "/networks/vtn/lab/SDMC_labscience/studies/CoVPN/CoVPN5001/assays/ADCC/misc_files/512-388-2021000250.txt"
    most_samples = pd.read_csv(most_samples_path, sep="\t")

    duke_samples_path = "/networks/vtn/lab/SDMC_labscience/studies/CoVPN/CoVPN5001/assays/ADCC/misc_files/512-388-2021000247.txt"
    duke_samples = pd.read_csv(duke_samples_path, sep="\t")

    most_samples.columns = [i.lower() for i in most_samples.columns]
    duke_samples.columns = [i.lower() for i in duke_samples.columns]

    most_samples_list = set(most_samples.global_id).difference(outputs.guspec)
    if len(most_samples_list)==0:
        print("All samples in 512-388-2021000250.txt are in outputs")

    not_accounted_for = set(outputs.guspec).difference(most_samples.global_id)
    if len(not_accounted_for) > 0:
        print(f"These samples aren't accounted for in that original manifest: {not_accounted_for}")

        duke_overlap = not_accounted_for.difference(duke_samples.global_id)
        if len(duke_overlap) == 0:
            print(f"The remaining samples not in the original manifest seem to be in this manifest sent to Duke (512-388-2021000247.txt). They must have shared them with the Ferrari Lab")


if __name__=="__main__":
    main()
