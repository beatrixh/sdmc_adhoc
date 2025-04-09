## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 04/09/2025
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
    input_data_path = '/trials/vaccine/p302/s001/qdata/LabData/AE_assays_pass-through/Anti-PEG/HVTN302 20487050 anti-Peg data 10March2025r1.xlsx'
    data = pd.read_excel(input_data_path, skiprows=2)

    # standardize column format
    data.columns = [i.lower().replace("\n", "").replace(" ", "_") for i in data.columns]

    # add guspec column
    data['guspec'] = data.sample_user_text_1.str[3:]
    data = data.rename(columns={'results_with_interpolated_titer':'result_as_submitted'})

    data['result_units'] = "Titer"
    data.loc[data.result_as_submitted.isin(['Negative Immunodepletion', 'Negative Screen']), 'result_units'] = "N/A"

    data['result_interpretation'] = "Positive"
    data.loc[data.result_as_submitted.isin(['Negative Immunodepletion', 'Negative Screen']), 'result_interpretation'] = "Negative"

    data['result'] = data['result_as_submitted']
    data.loc[data.result_as_submitted.isin(['Negative Titer (<50.000)']), 'result'] = "<50"

    datadir = '/trials/vaccine/p302/s001/qdata/LabData/AE_assays_pass-through/Anti-PEG/'
    original = datadir + '20487050 anti-Peg data 13May2024_with Titer_non QC_LB.xls'
    original = pd.read_excel(original)
    rerun_guspecs = original.loc[original['IR Result \nText']=='Positive Titer (>6400.000)']['Sample \nUser text 1'].str[3:].tolist()

    data['rerun'] = False
    data.loc[data.guspec.isin(rerun_guspecs), 'rerun'] = True

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

    # (pd.to_datetime(outputs.actual_sampling_date) != pd.to_datetime(outputs.drawdt)).sum()

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

    # merging on for consistency with original upload
    outputs['lab_internal_sample_receipt_date'] = '2024-02-22'
    outputs['analyte'] = "Anti-PEG_Antibodies"

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
        'result_as_submitted',
        'rerun',
        'assay_precision',
        'lab_internal_run_id',
        'lab_internal_sample_receipt_date',
        'sample_id_submitted',
        'sdmc_processing_datetime',
        'sdmc_data_receipt_datetime',
        'input_file_name'
    ]

    if len(list(set(reorder).symmetric_difference(outputs.columns))) > 0:
        print(set(reorder).symmetric_difference(outputs.columns))
        raise Exception("Oops; trying to drop or add columns")

    outputs = outputs[reorder]

    ## CHECKS ---------------------------------------------------------------------##

    ## everything with new run id previously was ">6400" --------------------------##
    previous = "/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN302/assays/AE_assays/Anti-PEG_total/misc_files/data_processing/HVTN302_Anti-PEG_total_CRL_interpolated_titers_processed_2025-04-04.txt"
    previous = pd.read_csv(previous, sep="\t")

    checkcols = ['guspec','lab_internal_run_id']
    check_runids = outputs[checkcols].merge(previous[checkcols], on='guspec', suffixes=('_new','_previous'))
    updated_runids = check_runids.loc[check_runids.lab_internal_run_id_new!=check_runids.lab_internal_run_id_previous].guspec.tolist()

    original = datadir + '20487050 anti-Peg data 13May2024_with Titer_non QC_LB.xls'
    original = pd.read_excel(original)
    rerun_guspecs = original.loc[original['IR Result \nText']=='Positive Titer (>6400.000)']['Sample \nUser text 1'].str[3:].tolist()

    assert(len(set(updated_runids).symmetric_difference(rerun_guspecs)) == 0), "updated run ids dont match >6400s"

    ## everything is contained in [old endpoint titer, old endpoint titer*2] -----##
    original['guspec'] = original['Sample \nUser text 1'].str[3:]
    original['result_numeric'] = original['IR Result \nText'].map({
        100: 100.,
        400: 400.,
        1600: 1600.,
        'Positive Titer (>6400.000)': 6400.,
        3200: 3200.,
        'Negative Screen': np.nan,
        'Negative Titer (<50.000)': 25.,
        800: 800.,
        200: 200.,
        'Negative Immunodepletion': np.nan,
        50: 50.
    })

    check_titers = outputs[['guspec','result_as_submitted']].merge(original[['guspec','result_numeric']], on='guspec', how='outer')
    check_titers.loc[check_titers.result_as_submitted.isin(['Negative Immunodepletion', 'Negative Screen']), 'result_as_submitted'] = np.nan
    check_titers.loc[check_titers.result_as_submitted=='Negative Titer (<50.000)', 'result_as_submitted'] = 25.
    check_titers.result_as_submitted = check_titers.result_as_submitted.astype(float)

    # this should only have gone down if it was rerun; eg if result numeric == '>6400'
    result_went_unexpectedly_down = check_titers.loc[(check_titers.result_as_submitted < check_titers.result_numeric) & (check_titers.result_numeric != 6400.)]
    assert len(result_went_unexpectedly_down) == 0, "result that we're rerun went down, oh no"

    # all results that weren't rerun should have gone up by no more than double
    result_went_unexpectedly_up_a_lot = check_titers.loc[(check_titers.result_as_submitted > check_titers.result_numeric *2) & (check_titers.result_numeric != 6400.)]
    assert len(result_went_unexpectedly_up_a_lot) == 0, "result that we're rerun went into a new endpoint titer, oh no"

    ## save outputs
    savedir = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN302/assays/AE_assays/Anti-PEG_total/misc_files/data_processing/'
    today = datetime.date.today().isoformat()
    fname = f"HVTN302_Anti-PEG_total_CRL_interpolated_titers_processed_{today}.txt"
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
