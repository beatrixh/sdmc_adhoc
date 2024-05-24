## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 03/08/2024
# Purpose:  Process/standardize RBM 302 IgE data
## ---------------------------------------------------------------------------##
import pandas as pd
import datetime
import os
import yaml

import sdmc_tools.process as sdmc
import sdmc_tools.constants as constants
# read in data ---------------------------------------------------------------##
def main():
    yaml_path = os.path.dirname(__file__) + "/paths.yaml"
    with open(yaml_path, 'r') as file:
        yaml_dict = yaml.safe_load(file)

    input_data_path = yaml_dict['input_data_path']

    # input_data_path = '/trials/vaccine/p302/s001/qdata/LabData/AE_assays_pass-through/IgE_RBM/HVTN302 Total IgE _ RBM_ 07MAr2024_130814 Results v1.xls'
    inputs = pd.read_excel(input_data_path, sheet_name="Custom Human Map", skiprows=0, header=0, index_col=None)

    # grab data vs metadata in chunks and reformat -------------------------------##

    # grab first chunk of metadata and reformat
    metadata1 = inputs.iloc[:5, [3, 7]]
    metadata1 = metadata1.T.rename(columns=metadata1.T.iloc[0]).iloc[1:].reset_index(drop=True)
    metadata1 = metadata1.drop(columns="Study ")

    # grab second chunk of metadata. this is all redundant so won't merge it on.
    metadata2 = inputs.iloc[[7],[3,7,8,11]]

    # grab third chunk of metadata
    metadata3 = inputs.iloc[9:15,[1,3]]
    metadata3 = metadata3.T.rename(columns=metadata3.T.iloc[0]).iloc[1:].reset_index(drop=True)

    # subset data out of input sheet
    data = inputs.iloc[17:296,[1,3]]

    # this should be 279
    print(f"len of data; should be 279: {len(data)}")

    # rename columns
    data.columns = ["sample_id_lab", "result"]

    # make sure all the guspecs start with fsq
    print(f"N guspecs starting with FSQ: {data.sample_id_lab.str.contains('FSQ').sum()}, len(data): {len(data)}")

    # create guspec column
    data["guspec"] = data.sample_id_lab.str[3:]

    # remove any whitespace
    data.guspec = data.guspec.replace(" ", "")

    # merge on lab-submitted metadata (we decided to drop #2)
    data = data.merge(metadata1, how="cross")
    data = data.merge(metadata3, how="cross")

    # standard processing ----------------------------------------------------##
    ldms = pd.read_csv(constants.LDMS_PATH_HVTN, usecols=constants.STANDARD_COLS)
    ldms = ldms.loc[ldms.lstudy==302.]
    ldms = ldms.loc[ldms.guspec.isin(data.guspec.unique().tolist())]

    metadata_dict = {
            'network': 'HVTN',
            'upload_lab_id': 'N4',
            'assay_lab_name': 'Rules-Based Medicine (RBM)',
            'instrument': 'Luminex',
            'assay_type': 'DiscoveryMAP',
            'specrole': 'Sample',
    }

    outputs = sdmc.standard_processing(
        input_data=data,
        input_data_path=input_data_path,
        guspec_col="guspec",
        network="hvtn",
        metadata_dict=metadata_dict,
        ldms=ldms
    )

    # cleanup ----------------------------------------------------------------##
    outputs = outputs.rename(columns = {
        'analytes': 'analyte',
        'rbm_ldd': 'llod',
        'rbm_lloq': 'lloq',
        'rbm_serum_low_range': 'reference_range_low',
        'rbm_serum_high_range': 'reference_range_high',
    })

    # got rid of this column; it was all NaN
    outputs = outputs.drop(columns="group_designation")

    reorder = ['network', 'protocol', 'specrole', 'guspec', 'ptid', 'visitno',
               'drawdt', 'spectype', 'spec_primary', 'spec_additive',
               'spec_derivative', 'upload_lab_id', 'assay_lab_name', 'assay_type',
               'instrument', 'analyte', 'units', 'result', 'llod', 'lloq',
               'reference_range_low', 'reference_range_high', 'sample_id_lab',
               'rbm_order', 'customer', 'investigator',
               'sdmc_processing_datetime', 'sdmc_data_receipt_datetime',
               'input_file_name']

    col_diff = set(reorder).symmetric_difference(outputs.columns)
    if len(list(col_diff)) > 0:
        raise Exception(f"The following cols being dropped or added: {col_diff}")

    # reorder columns
    outputs = outputs[reorder]

    # save -------------------------------------------------------------------##
    today = datetime.datetime.today().strftime('%Y-%m-%d')
    save_fname = "DRAFT_" + yaml_dict["output_prefix"] + f"_{today}.txt"

    outputs.to_csv(yaml_dict["savedir"] + save_fname, sep="\t", index=False)

if __name__ == '__main__':
    main()
