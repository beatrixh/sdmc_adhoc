## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 07/16/2024
# Purpose: Process ICABA data from Ferrari lab
## ---------------------------------------------------------------------------##
import pandas as pd
import datetime, os
import yaml
import sdmc_tools.process as sdmc
import sdmc_tools.constants as constants

def main():
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

    ldms = pd.read_csv(
        constants.LDMS_PATH_HVTN,
        usecols=constants.STANDARD_COLS,
        dtype=constants.LDMS_DTYPE_MAP
        )

    if len(set(input_data.guspec).difference(ldms.guspec)) > 0:
        raise Exception("input data has guspecs missing from ldms")

    ldms = ldms.loc[ldms.guspec.isin(input_data.guspec.unique())]

    outputs = sdmc.standard_processing(
        input_data=input_data,
        input_data_path=input_data_path,
        guspec_col='guspec',
        network='hvtn',
        metadata_dict=metadata_dict,
        ldms=ldms
    )

    # fix filepath sources for data we receieved separately
    addendum_guspecs = addendum_data.guspec.str[3:].unique().tolist()
    outputs.loc[outputs.guspec.isin(addendum_guspecs), 'sdmc_data_receipt_datetime'] = datetime.datetime.fromtimestamp(os.path.getmtime(addendum_data_path)).replace(microsecond=0).isoformat()
    outputs.loc[outputs.guspec.isin(addendum_guspecs), 'input_file_name'] = addendum_data_path.split("/")[-1]

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
        'log_titer',
        'titer_calculated_by_sdmc',
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

if __name__=="__main__":
    main()
