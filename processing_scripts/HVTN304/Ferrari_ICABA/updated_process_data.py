## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 05/14/2024
# Purpose:  - Process ICABA data from Ferrari lab
## ---------------------------------------------------------------------------##
import pandas as pd
import numpy as np
import datetime
import yaml
import os

import sdmc_tools.refactor_process as sdmc
import sdmc_tools.constants as constants
## ---------------------------------------------------------------------------##

def main():
    with open('paths.yaml', 'r') as file:
        yaml_dict = yaml.safe_load(file)

    input_data_path = yaml_dict['input_data_path'][0]
    metadata_path = yaml_dict['input_data_path'][1]

    update_data_path = "/trials/vaccine/p304/s001/qdata/LabData/ICABA_pass-through/Ferrari_HVTN 304_ICABA_Analysis_23JUL2024.csv"
    update_metadata_path = "/trials/vaccine/p304/s001/qdata/LabData/ICABA_pass-through/Ferrari_HVTN 304_ICABA_Metadata_23JUL2024.csv"

    ldms = pd.read_csv(
        constants.LDMS_PATH_HVTN,
        usecols=constants.STANDARD_COLS,
        dtype=constants.LDMS_DTYPE_MAP
    )
    ldms = ldms.loc[ldms.lstudy==304.]

    def prep_data(input_data_path, metadata_path):
        data = pd.read_csv(input_data_path)
        metadata = pd.read_csv(metadata_path)

        data = data.dropna(how='all', axis=1)

        data = data.rename(columns={"sample_id":"guspec",
                                "visit":"visitno",
                                'Mock%IgG+': 'result_Mock%IgG+',
                                '%p24+IgG+': 'result_%p24+IgG+',
                                'MockSubtracted %p24+IgG+': 'result_MockSubtracted %p24+IgG+'
                               })
        metadata = metadata.rename(columns={'sample_id':'guspec',
                                        'visit':'visitno'})
        data = data.merge(metadata,
                          on=["guspec", "ptid", "visitno", "dilution"],
                          how='outer'
        )

        data = data.drop(columns=["visitno", "ptid", "analysis_file_name"])
        data = data.rename(columns={'xml_file_name':'lab_xml_file_name'})

        md = {
            'network': 'HVTN',
            'upload_lab_id': 'GF',
            'assay_lab_name': 'Ferrari Lab',
            'instrument': 'Fortessa Flow Cytometer',
            'assay_type': 'ICABA',
            'specrole': 'Sample',
            'result_units': 'Percent',
        }

        outputs = sdmc.standard_processing(
            input_data=data,
            input_data_path=input_data_path,
            additional_input_paths={'input_metadata_file_name':metadata_path},
            guspec_col="guspec",
            network="HVTN",
            metadata_dict=md,
            ldms=ldms
        )

        return outputs

    original_outputs = prep_data(input_data_path, metadata_path)
    new_outputs = prep_data(update_data_path, update_metadata_path)

    # replace data with new data
    ptid_to_replace = new_outputs.ptid.unique().tolist()
    original_outputs = original_outputs.loc[~original_outputs.ptid.isin(ptid_to_replace)]

    outputs = pd.concat([original_outputs, new_outputs])

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
        'result_%p24+igg+',
        'result_mock%igg+',
        'result_mocksubtracted_%p24+igg+',
        'result_units',
        'virus_id',
        'virus_stock',
        'imc_prep',
        'assay_date',
        'fcs_file_name',
        'isotype',
        'lab_xml_file_name',
        'operator',
        'sdmc_processing_datetime',
        'sdmc_data_receipt_datetime',
        'input_file_name',
        'input_metadata_file_name'
    ]

    mismatch = set(outputs.columns).symmetric_difference(reorder)
    if len(mismatch) > 0:
        raise Exception(f"trying to drop or add columns on reorder: {mismatch}")

    outputs = outputs[reorder]
    outputs.columns = [i.lower().replace(" ","_") for i in outputs.columns]
    outputs['input_file_name'] = f"{input_data_path.split('/')[-1]}, {metadata_path.split('/')[-1]}"

    today = datetime.date.today().isoformat()
    savepath = yaml_dict["savedir"] + "DRAFT_" + yaml_dict["output_prefix"] + f"_{today}.txt"
    outputs.to_csv(savepath, sep="\t", index=False)

if __name__ == '__main__':
    main()
