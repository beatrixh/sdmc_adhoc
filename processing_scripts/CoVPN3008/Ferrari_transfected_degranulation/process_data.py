## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 06/10/2024
# Purpose:  - Process Tranfected CD107a Degranulation data from Ferrari lab
## ---------------------------------------------------------------------------##
import pandas as pd
import datetime, os
import yaml
import sdmc_tools.process as sdmc
import sdmc_tools.constants as constants

def main():

    # pull in all data and merge
    dirpath = "/trials/covpn/p3008/s001/qdata/LabData/Transfected_CD107a_degranulation_pass-through/"
    ba2_analysis_fname = "Ferrari_CoVPN3008_BA2.86_Transfected Cell Degranulation_Analysis_24JAN2024.csv"
    ba2_metadata_fname = "Ferrari_CoVPN3008_BA2.86_Transfected Cell Degranulation_Metadata_24JAN2024.csv"

    ba4_analysis_fname = "Ferrari_CoVPN3008_BA4_5_Transfected Cell Degranulation_Analysis_24JAN2024.csv"
    ba4_metadata_fname = "Ferrari_CoVPN3008_BA4_5_Transfected Cell Degranulation_Metadata_24JAN2024.csv"

    g6_analysis_fname = "Ferrari_CoVPN3008_G614_Transfected Cell Degranulation_Analysis_24JAN2024.csv"
    g6_metadata_fname = "Ferrari_CoVPN3008_G614_Transfected Cell Degranulation_Metadata_24JAN2024.csv"

    ba2 = pd.read_csv(dirpath + ba2_analysis_fname)
    ba2_meta = pd.read_csv(dirpath + ba2_metadata_fname)
    ba2 = ba2.merge(ba2_meta, on=["sample_id", "ptid", "visit", "dilution"], how="outer")

    ba4 = pd.read_csv(dirpath + ba4_analysis_fname)
    ba4_meta = pd.read_csv(dirpath + ba4_metadata_fname)
    ba4 = ba4.merge(ba4_meta, on=["sample_id", "ptid", "visit", "dilution"], how="outer")

    g6 = pd.read_csv(dirpath + g6_analysis_fname)
    g6_meta = pd.read_csv(dirpath + g6_metadata_fname)
    g6 = g6.merge(g6_meta, on=["sample_id", "ptid", "visit", "dilution"], how="outer")

    ba2.columns = [i.lower().replace(" ", "_").replace("%", "pct_") for i in ba2.columns]
    ba4.columns = [i.lower().replace(" ", "_").replace("%", "pct_") for i in ba4.columns]
    g6.columns = [i.lower().replace(" ", "_").replace("%", "pct_") for i in g6.columns]

    ba2['antigen'] = 'BA2.86'
    ba4['antigen'] = 'BA4_5'
    g6['antigen'] = 'G614'

    data = pd.concat([ba2, ba4, g6])
    data = data.rename(columns={
        'sample_id':'guspec',
        'ptid':'submitted_ptid',
        'visit':'submitted_visit',
        'infected/transfected': 'antigen_value'
    })
    data["antigen_value"] = "Transfected"

    # read in ldms
    ldms = pd.read_csv(
        constants.LDMS_PATH_COVPN,
        usecols=constants.STANDARD_COLS,
        dtype=constants.LDMS_DTYPE_MAP
    )
    ldms = ldms.loc[ldms.lstudy==3008.]
    ldms.txtpid = ldms.txtpid.astype(str)
    ldms = ldms.loc[ldms.guspec.isin(ba2.sample_id.unique())].drop_duplicates()

    # metadata from sdmc background materials
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

    # standard processing
    outputs = sdmc.standard_processing(
        input_data=data,
        input_data_path=dirpath + ba4_analysis_fname,
        guspec_col='guspec',
        network='covpn',
        metadata_dict=metadata_dict,
        ldms=ldms
    )

    # input file correction -- sdmc.standard_processing cant handle multiple inputs
    input_file_map = {
        'BA2.86': f'{ba2_analysis_fname}, {ba2_metadata_fname}',
        'BA4_5': f'{ba4_analysis_fname}, {ba4_metadata_fname}',
        'G614': f'{g6_analysis_fname}, {g6_metadata_fname}'
    }
    outputs['input_file_name'] = outputs.antigen.map(input_file_map)
    outputs = outputs.drop(columns=["submitted_ptid", "submitted_visit", "infected/transfected"])

    # reorder outputs
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
        'antigen',
        'dilution',
        'infected_pct_cd107a+_nk_cells',
        'mock_pct_cd107a+_nk_cells',
        'background_subtracted_pct_cd107a+_nk_cells',
        'result_units',
        'assay_date',
        'operator',
        'analysis_file_name',
        'fcs_file_name',
        'xml_file_name',
        'sdmc_processing_datetime',
        'sdmc_data_receipt_datetime',
        'input_file_name',
    ]

    # set(reorder).symmetric_difference(outputs.columns)
    outputs = outputs[reorder]

    # save outputs
    savedir = '/networks/vtn/lab/SDMC_labscience/studies/CoVPN/CoVPN3008/assays/tranfected_cell_cd107a_degranulation/misc_files/data_processing/'
    today = datetime.datetime.today().strftime('%Y-%m-%d')
    fname = f"CoVPN3008_Ferrari_ADCC_Degranulation_processed_{today}.txt"

    outputs.to_csv(savedir + "DRAFT_" + fname, index=False, sep="\t")

if __name__ == '__main__':
    main()
