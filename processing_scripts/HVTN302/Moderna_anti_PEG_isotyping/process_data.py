## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 10-29-2024
# Purpose: Process Anti-PEG isotyping (adding on data from anti-peg total negatives)
## ---------------------------------------------------------------------------##
import pandas as pd
import numpy as np
import datetime, os
import yaml
import sdmc_tools.process as sdmc
import sdmc_tools.constants as constants

def main():
    ## read in ldms ----------------------------------------------------------##
    ldms = pd.read_csv(
        # constants.LDMS_PATH_HVTN,
        '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN302/specimens/ldms_feed/hvtn.ldms302.20241030.csv',
        usecols=constants.STANDARD_COLS,
        dtype=constants.LDMS_DTYPE_MAP
    )
    ldms = ldms.loc[ldms.lstudy==302.]

    ## -----------------------------------------------------------------------##
    input_data_path = "/trials/vaccine/p302/s001/qdata/LabData/AE_assays_pass-through/Anti-PEG/BMA_2024_0031 Final Data.xlsx"
    input_data = pd.read_excel(input_data_path)

    addendum_data_path = "/trials/vaccine/p302/s001/qdata/LabData/AE_assays_pass-through/Anti-PEG/BMA_2024_0031 Final Data_Addendum_12Jul24.xlsx"
    addendum_data = pd.read_excel(addendum_data_path)

    input_data = pd.concat([input_data, addendum_data])
    input_data.columns = [i.lower().replace(" ","_") for i in input_data.columns]
    input_data.guspec = input_data.guspec.str[3:]

    input_data = input_data.rename(columns={
        'specrole':'sample_id_lab',
        'watson_run_id':'run_id'
    })
    def exponentiate(logtiter):
        try:
            return 10**float(logtiter)
        except:
            return logtiter
    input_data["titer_calculated_by_sdmc"] = input_data.log_titer.apply(lambda x: exponentiate(x))

    # add dilution ranges
    dilution_path = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN302/assays/AE_assays/Anti-PEG_isotyping/misc_files/dilution_ranges.xlsx'
    dilutions = pd.read_excel(dilution_path)
    dilutions = dilutions.set_index("Dilution")

    mins = {ig: dilutions.loc['Dil 1', ig] for ig in ['IgG1', 'IgG3', 'IgG4', 'IgM', 'IgE']}
    maxs = {ig: dilutions.loc['Dil 8', ig] for ig in ['IgG1', 'IgG3', 'IgG4', 'IgM', 'IgE']}

    input_data["dilution_range_min"] = input_data.isotype.map(mins)
    input_data["dilution_range_max"] = input_data.isotype.map(maxs)

    metadata_dict = {'network': 'HVTN',
                     'upload_lab_id': 'N4',
                     'assay_lab_name': 'Moderna',
                     'instrument': 'MSD',
                     'lab_software_version': 'MSD Discovery Workbench, Excel, JMP',
                     'assay_type': 'Anti-PEG Isotyping',
                     'specrole': 'Sample',
                     'assay_precision': 'Semi-Quantitative'}

    outputs = sdmc.standard_processing(
        input_data=input_data,
        input_data_path=input_data_path,
        guspec_col='guspec',
        network='hvtn',
        metadata_dict=metadata_dict,
        ldms=ldms.loc[ldms.guspec.isin(input_data.guspec)]
    )

    # fix filepath sources for data we receieved separately
    addendum_guspecs = addendum_data.guspec.str[3:].unique().tolist()
    outputs.loc[outputs.guspec.isin(addendum_guspecs), 'sdmc_data_receipt_datetime'] = datetime.datetime.fromtimestamp(os.path.getmtime(addendum_data_path)).replace(microsecond=0).isoformat()
    outputs.loc[outputs.guspec.isin(addendum_guspecs), 'input_file_name'] = addendum_data_path.split("/")[-1]

    # add 'negative type' to differentiate from those screened out in anti-peg total
    outputs['negative_type'] = "N/A"
    outputs.loc[outputs.log_titer=="No Titer",'negative_type'] = "No Titer"
    outputs.loc[outputs.log_titer=="No Titer",'titer_calculated_by_sdmc'] = outputs.loc[outputs.log_titer=="No Titer",'dilution_range_min'] / 2
    outputs.titer_calculated_by_sdmc = outputs.titer_calculated_by_sdmc.astype(float)
    outputs.loc[outputs.log_titer=="No Titer",'log_titer'] = np.log10(outputs.loc[outputs.log_titer=="No Titer",'titer_calculated_by_sdmc'])

    ## -----------------------------------------------------------------------##
    total_path = "/trials/vaccine/p302/s001/qdata/LabData/AE_assays_pass-through/Anti-PEG/processed_by_sdmc/HVTN302_Anti-PEG_total_CRL_processed_2024-07-02.txt"
    total = pd.read_csv(total_path, sep="\t")
    negatives = total.loc[total.result_interpretation=="Negative"]

    negatives = negatives[['guspec','sample_id_submitted','result']].drop_duplicates()

    assay_metadata = input_data[['isotype','cutpoint','dilution_range_min','dilution_range_max','sensitivity']].drop_duplicates()
    negatives = negatives.merge(assay_metadata, how='cross')
    negatives = negatives.rename(columns={
        'result':'negative_type',
        'sample_id_submitted': 'sample_id_lab'
    })

    metadata_dict_negatives = {
        'network': 'HVTN',
        'upload_lab_id': 'N4',
        'assay_lab_name': 'Moderna',
        'instrument': 'N/A (value estimated)',
        'lab_software_version': 'N/A (value estimated)',
        'assay_type': 'Anti-PEG Isotyping',
        'specrole': 'Sample',
        'assay_precision': 'Semi-Quantitative',
        'run_id': 'N/A (merged from anti-PEG total)'
    }

    outputs_negatives = sdmc.standard_processing(
        input_data=negatives,
        input_data_path=total_path,
        guspec_col='guspec',
        network='hvtn',
        metadata_dict=metadata_dict_negatives,
        ldms=ldms.loc[ldms.guspec.isin(negatives.guspec)]
    )

    outputs_negatives['titer_calculated_by_sdmc'] = outputs_negatives.dilution_range_min / 2
    outputs_negatives['log_titer'] = np.log10(outputs_negatives.titer_calculated_by_sdmc)

    ## -----------------------------------------------------------------------##
    col_diff = set(outputs.columns).symmetric_difference(outputs_negatives.columns)
    if len(col_diff) > 0:
        raise Exception("Different columns in outputs and outputs_negatives!")

    outputs = pd.concat([outputs, outputs_negatives])
    outputs.log_titer = outputs.log_titer.astype(float)

    log_check = np.abs(10**outputs.log_titer - outputs.titer_calculated_by_sdmc).max()
    if log_check > 0:
        raise Exception("log_titer isn't the log of titer! go double check")

    # these should cancel out
    est_titer = np.abs(outputs.loc[outputs.negative_type!="N/A"].dilution_range_min - outputs.loc[outputs.negative_type!="N/A"].titer_calculated_by_sdmc*2).sum()
    if est_titer > 0:
        raise Exception("half of min dilution calc for negatives done incorrectly!")

    # all of the negatives should come from the 49 that were screened out in anti-peg total, plus the negative titers from isotyping
    check_count = len(outputs.loc[outputs.negative_type!="N/A"]) - (49*5 + 713)
    if check_count != 0:
        raise Exception("The count of negatives doesn't match expected!")

    log_titer_source_map = {
        'No Titer': 'Estimated as log of half the min dilution per isotype',
        'N/A': 'Lab-calculated using logistic fit to log dilution',
        'Negative Screen': 'Estimated as log of half the min dilution per isotype',
        'Negative Immunodepletion': 'Estimated as log of half the min dilution per isotype',
    }

    outputs['log_titer_source'] = outputs.negative_type.map(log_titer_source_map)

    ## -----------------------------------------------------------------------##
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
        'lab_software_version',
        'isotype',
        'dilution_range_min',
        'dilution_range_max',
        'cutpoint',
        'sensitivity',
        'negative_type',
        'log_titer',
        'titer_calculated_by_sdmc',
        'log_titer_source',
        'assay_precision',
        'sample_id_lab',
        'run_id',
        'sdmc_processing_datetime',
        'sdmc_data_receipt_datetime',
        'input_file_name'
    ]

    mismatch = set(reorder).symmetric_difference(outputs.columns)
    if len(mismatch) > 0:
        raise Exception(f"trying to drop or add columns: {mismatch}")

    outputs = outputs[reorder]

    today = datetime.datetime.today().strftime('%Y-%m-%d')
    savedir = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN302/assays/AE_assays/Anti-PEG_isotyping/misc_files/data_processing/'
    fname = f"DRAFT_HVTN302_Anti-PEG_isotyping_Moderna_processed_{today}.txt"
    outputs.to_csv(savedir + fname, index=False, sep="\t")

    ptid_visitno_isotype = pd.pivot_table(outputs, index=['ptid','visitno'], columns='isotype', aggfunc='count', fill_value=0)['titer_calculated_by_sdmc']
    ptid_visitno_isotype.to_excel(savedir + "HVTN302_antipeg_isotyping_")


if __name__=="__main__":
    main()
