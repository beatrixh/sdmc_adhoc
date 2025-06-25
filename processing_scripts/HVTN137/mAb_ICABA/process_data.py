## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 06/25/2025
# Purpose:  Ad hoc data processing of mAb ICABA data from Ferrari lab
## ---------------------------------------------------------------------------##
import pandas as pd
import numpy as np
import datetime, os
import yaml
import sdmc_tools.process as sdmc
import sdmc_tools.constants as constants
import sdmc_tools.access_ldms as access_ldms

def main():
    ## read in data + metadata -----------------------------------------------##
    analysis_path = '/trials/vaccine/p137/s001/qdata/LabData/Ferrari_ICABA_pass-through/Ferrari_HVTN 137_mAbs_ICABA_Analysis_20JUN25.csv'
    data = pd.read_csv(analysis_path)

    metadata_path = '/trials/vaccine/p137/s001/qdata/LabData/Ferrari_ICABA_pass-through/Ferrari_HVTN 137_mAbs_ICABA_Metadata_20JUN25.csv'
    metadata = pd.read_csv(metadata_path)

    ## rename columns and merge ----------------------------------------------##
    data = data.rename(columns={"sample_id":"mab_id",
                            'IMC':'virus_id',
                            'Mock%IgG+': 'result_Mock%IgG+',
                            '%p24+IgG+': 'result_%p24+IgG+',
                            'MockSubtracted %p24+IgG+': 'result_MockSubtracted %p24+IgG+',
                            '%GFP+IgG+':'result_%GFP+IgG+',
                            'MockSubtracted %GFP+IgG+':'result_MockSubtracted %GFP+IgG+'
                           })

    metadata = metadata.rename(columns={
        'sample_id':'mab_id',
        'xml_file_name':'lab_xml_file_name'
    })
    metadata = metadata.dropna(how='all')

    for df in [metadata, data]:
        df.columns = [i.lower() for i in df.columns]

    metadata.concentration_mg = metadata.concentration_mg.astype(int)
    data = data.merge(metadata, on = ['mab_id','virus_id','concentration_mg'], how = 'outer')

    data = data.drop(columns=["analysis_file_name","unnamed: 8"])

    data = data.rename(columns={
        'concentration_mg':'mab_concentration'
    })
    data['mab_concentration_units'] = 'mg'

    ## merge on hand-pulled metadata + std formatting ------------------------##
    md = {
        'network': 'HVTN',
        'protocol': 137,
        'upload_lab_id': 'GF',
        'assay_lab_name': 'Ferrari Lab',
        'instrument': 'Fortessa Flow Cytometer',
        'assay_type': 'ICABA',
        'specrole': 'Sample',
        'spectype': 'Monoclonal antibody',
        'result_units': 'Percent',
    }

    outputs = sdmc.processing_minus_ldms(
        input_data=data,
        metadata_dict=md,
        input_data_path=analysis_path,
        additional_input_paths={'input_metadata_file_name':metadata_path},
        cols_to_lower=True
    )

    reorder = [
        'network',
        'protocol',
        'specrole',
        'spectype',
        'mab_id',
        'upload_lab_id',
        'assay_lab_name',
        'assay_type',
        'instrument',
        'assay_date',
        'isotype',
        'mab_concentration',
        'mab_concentration_units',
        'result_%gfp+igg+',
        'result_%p24+igg+',
        'result_mock%igg+',
        'result_mocksubtracted_%gfp+igg+',
        'result_mocksubtracted_%p24+igg+',
        'result_units',
        'virus_id',
        'virus_stock',
        'imc_prep',
        'operator',
        'lab_xml_file_name',
        'fcs_file_name',
        'sdmc_processing_datetime',
        'sdmc_data_receipt_datetime',
        'input_file_name',
        'input_metadata_file_name',
    ]

    assert len(set(reorder).symmetric_difference(outputs.columns)) == 0
    outputs = outputs[reorder]

    ## checks ----------------------------------------------------------------##
    check_p24 = outputs['result_%p24+igg+'] - outputs['result_mock%igg+']
    check_p24[check_p24 < 0] = 0

    if np.abs(check_p24 - outputs['result_mocksubtracted_%p24+igg+']).max() > 1e-10:
        raise Exception("Subtraction error")

    assert len(outputs.loc[(
        (outputs['result_%gfp+igg+'].isna() & (outputs['result_mocksubtracted_%gfp+igg+'].notna())) |
        (outputs['result_%gfp+igg+'].notna() & (outputs['result_mocksubtracted_%gfp+igg+'].isna()))
    )]) == 0

    check_gfp = outputs['result_%gfp+igg+'] - outputs['result_mock%igg+']
    check_gfp[check_gfp < 0] = 0

    if np.abs(check_gfp - outputs['result_mocksubtracted_%gfp+igg+']).max() > 1e-10:
        raise Exception("Subtraction error")

    ## save to txt file ------------------------------------------------------##
    savedir = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN137/assays/ICABA/misc_files/data_processing_mAb/'
    today = datetime.date.today().isoformat()
    outputs.to_csv(savedir + f"HVTN137_mAb_ICABA_processed_{today}.txt", sep="\t", index=False)

    ## pivot summary ---------------------------------------------------------##
    pivot_summary = pd.pivot_table(
        outputs,
        index=['mab_id','mab_concentration'],
        columns=['virus_id'],
        aggfunc='count'
    )[['result_mocksubtracted_%p24+igg+']]

    pivot_summary.to_excel(savedir + "HVTN137_mAb_ICABA_pivot_summary.xlsx")

if __name__=="__main__":
    main()
