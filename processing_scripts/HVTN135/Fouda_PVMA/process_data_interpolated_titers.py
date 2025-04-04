## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 04/04/2025
# Purpose:  - Process updated CRL total anti-PEG data processing
## ---------------------------------------------------------------------------##
import pandas as pd
import numpy as np
import datetime, os
import yaml
import sdmc_tools.process as sdmc
import sdmc_tools.constants as constants

## Read in data --------------------------------------------------------------##
def main():
    input_data_path = '/trials/vaccine/p302/s001/qdata/LabData/AE_assays_pass-through/Anti-PEG/HVTN302 20487050 anti-Peg data 10March2025.xls'
    data = pd.read_excel(input_data_path, skiprows=2)

    # standardize column format
    data.columns = [i.lower().replace("\n", "").replace(" ", "_") for i in data.columns]

    # add guspec column
    data['guspec'] = data.sample_user_text_1.str[3:]
    data = data.rename(columns={'results_with_extrapolated_titer':'result_as_submitted'})

    data['result_units'] = "Titer"
    data.loc[data.result_as_submitted.isin(['Negative Immunodepletion', 'Negative Screen']), 'result_units'] = "N/A"

    data['result_interpretation'] = "Positive"
    data.loc[data.result_as_submitted.isin(['Negative Immunodepletion', 'Negative Screen']), 'result_interpretation'] = "Negative"

    ## ldms
    ldms_dir = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN302/specimens/ldms_feed/'
    ldms = pd.read_csv(ldms_dir + np.sort(os.listdir(ldms_dir))[-1], usecols=constants.STANDARD_COLS, dtype=constants.LDMS_DTYPE_MAP)
    ldms = ldms.loc[ldms.guspec.isin(data.guspec)]

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

    (pd.to_datetime(outputs.actual_sampling_date) != pd.to_datetime(outputs.drawdt)).sum()

    ## checks ##
    if sum(outputs.subject.astype(int)!=outputs.ptid.astype(int)) > 0:
        raise Exception("ptids don't match ldms")

    ## drop columns
    outputs = outputs.drop(columns=[
        'actual_sampling_date',
        'actual_sampling_date/time',
        'subject',
        'study_title',
        'custom_id',
        'day_nominal',
        'nominal_time',
        'nominal_time_(design)',
        'sample_name',
        'study_day'
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
        'assay_precision',
        'result_as_submitted',
        'result_interpretation',
        'result_units',
        'lab_internal_run_id',
        'sample_id_submitted',
        'sdmc_processing_datetime',
        'sdmc_data_receipt_datetime',
        'input_file_name'
    ]

    if len(list(set(reorder).symmetric_difference(outputs.columns))) > 0:
        print(set(reorder).symmetric_difference(outputs.columns))
        raise Exception("Oops; trying to drop or add columns")

    outputs = outputs[reorder]

    ## save outputs
    savedir = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN302/assays/AE_assays/Anti-PEG_total/misc_files/data_processing/'
    today = datetime.date.today().isoformat()
    fname = f"VTN302_Anti-PEG_total_CRL_interpolated_titers_processed_{today}.txt"
    outputs.to_csv(savedir + fname, sep="\t", index=False)

    ## check for completeness against manifest
    manifest = pd.read_csv(savedir + "shipping_manifest.txt", sep="\t")
    outputs['core'] = outputs.guspec.str[:-4]
    manifest['core'] = manifest.GLOBAL_ID.str[:-4]

    # we have an exact match for guspec cores
    core_mismatch = set(manifest.core).symmetric_difference(outputs.core)
    print(f"exact match for guspec cores. core_mismatch: {core_mismatch}")

    # all outputs from the lab are in the manifest
    missing_from_manifest = set(outputs.guspec).difference(manifest.GLOBAL_ID)
    print(f"all samples from lab are in the manifest")
    print(f"sample from lab missing from manifest: {missing_from_manifest}")

    # there are exactly three aliquots per guspec core in the manifest
    aliquots_per_core = manifest.groupby('core')['GLOBAL_ID'].nunique().value_counts()
    print(f"all cores have three aliquots: {aliquots_per_core}")

    pd.pivot_table(data=outputs, values='result', index='ptid', columns='visitno', aggfunc='count', fill_value=0).to_csv(savedir + "ptid_visitno_summary.csv")
    outputs.guspec.value_counts().to_frame().to_csv(savedir + "guspec_summary.csv")

def pivot():
    pd.pivot_table(data=outputs, values='result', index='ptid', columns='visitno', aggfunc='count', fill_value=0).to_csv(savedir + "ptid_visitno_summary.csv")
    outputs.guspec.value_counts().to_frame().to_csv(savedir + "guspec_summary.csv")

if __name__ == '__main__':
    main()
