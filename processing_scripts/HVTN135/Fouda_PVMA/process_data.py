## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 03/31/2025
# Purpose:  Ad hoc data processing of PVMA data from Fouda Lab
## ---------------------------------------------------------------------------##
import pandas as pd
import numpy as np
import datetime
import yaml
import os

import sdmc_tools.process as sdmc
import sdmc_tools.constants as constants
## ---------------------------------------------------------------------------##

def main():
    # compare the previous upload to this one --------------------------------##
    # input_data_path = '/trials/vaccine/p135/s001/qdata/LabData/PVMA_pass-through/20250321_HVTN135_DataSummary_V3_MFI and Concentration.xlsx'
    # data = pd.read_excel(input_data_path)

    # old_data_path = '/trials/vaccine/p135/s001/qdata/LabData/PVMA_pass-through/20240822_HVTN135_DataSummary_V2.xlsx'
    # old = pd.read_excel(old_data_path)

    # # verify no changes from prev file
    # set(data.columns).difference(old.columns)
    # set(old.columns).difference(data.columns)
    # data[['Data Type','MFI < 100']].drop_duplicates()

    # looks like the lab-internal labels don't match from 33 onward (off-by-1 error)

    # input_data_path = '/trials/vaccine/p135/s001/qdata/LabData/PVMA_pass-through/20250321_HVTN135_DataSummary_V3_MFI and Concentration.xlsx'
    # data_mfi = pd.read_excel(input_data_path).rename(columns=rename)
    # data_iu = pd.read_excel(input_data_path, sheet_name="Analyzed Concentrations").rename(columns=rename)

    # data_mfi['id'] = data_mfi.lab_internal_sample_id.str.split(" ", expand=True)[1].astype(int)
    # data_iu['id'] = data_iu.lab_internal_sample_id.str.split(" ", expand=True)[1].astype(int)

    # test = data_mfi[['guspec','id']].merge(data_iu[['guspec','id']], on='guspec', how='outer')
    # test.loc[test.id_x!=test.id_y].sort_values(by='id_x')

    # read in and cast MFI data to long --------------------------------------##
    input_data_path = '/trials/vaccine/p135/s001/qdata/LabData/PVMA_pass-through/20250321_HVTN135_DataSummary_V3_MFI and Concentration.xlsx'
    data_mfi = pd.read_excel(input_data_path)

    # cast MFI data to long
    rename = {
        'GLOBAL_ID':'guspec',
        'Tube Label (WCM Internal Only)': 'lab_internal_sample_id',
        'Plate #': 'plate_number',
        'internal_qc':'lab_internal_qc'
    }
    data_mfi = data_mfi.rename(columns=rename)

    result_cols = [i for i in data_mfi.columns if '(' in i]
    non_result_cols = list(set(data_mfi.columns).difference(result_cols))
    data_mfi = data_mfi.melt(id_vars = non_result_cols, value_vars = result_cols, var_name='antigen_string', value_name='result_mfi')

    mfi = data_mfi.antigen_string.str.rpartition(" ", expand=True)
    mfi.columns = ['antigen', 'na', 'beadset_num']
    data_mfi = pd.concat([data_mfi, mfi], axis=1)
    data_mfi.beadset_num = data_mfi.beadset_num.str[1:-1].astype(int)

    # read in and cast IU data to long ---------------------------------------##
    data_iu = pd.read_excel(input_data_path, sheet_name="Analyzed Concentrations")
    data_iu = data_iu.rename(columns=rename)

    result_cols = [i for i in data_iu.columns if '(' in i]
    non_result_cols = list(set(data_iu.columns).difference(result_cols))
    data_iu = data_iu.melt(id_vars = non_result_cols, value_vars = result_cols, var_name='antigen_string', value_name='result_concentration')

    data_iu['MFI < 100'] = np.nan
    data_iu = data_iu.drop(columns='Blank')

    concentrations = data_iu.antigen_string.str.split("(", expand=True)
    concentrations.columns = ['antigen', 'result_concentration_units']
    concentrations.antigen = concentrations.antigen.str.strip()
    concentrations.result_concentration_units = concentrations.result_concentration_units.str[:-1].str.strip()

    data_iu = pd.concat([data_iu, concentrations], axis=1)

    # verify these match
    # set(data_iu.guspec).symmetric_difference(data_mfi.guspec)
    # data_iu[['guspec','antigen']].drop_duplicates().shape[0] == data_iu.shape[0]
    # data_mfi[['guspec','antigen']].drop_duplicates().shape[0] == data_mfi.shape[0]

    # merge rlu and concentrations together ----------------------------------##
    mergecols = [
        'guspec',
        'antigen',
        'plate_number',
        'DER',
        'Assay Tech',
        'Assay Lab',
        'Assay Date',
        'Internal QC',
        'MFI < 100',
        'Sample Dilution'
    ]

    dropcols = ['antigen_string','Data Type', 'Upload Date']
    data = data_mfi.drop(columns=dropcols).merge(data_iu.drop(columns=dropcols),
                                                         on = mergecols,
                                                         how='outer',
                                                         suffixes=('_20240822_upload','_20250321_upload')
                                                        )

    data = data.drop(columns='na')
    antigen_to_concentration_units = data[['antigen','result_concentration_units']].drop_duplicates()

    # read in std curve data then merge on -----------------------------------##

    std_curve = pd.read_excel(input_data_path, sheet_name="Antigen Standards")

    # consolidate header rows
    std_curve.columns = [f"{i} {j}" if 'Unnamed' not in i else j for (i,j) in zip(std_curve.columns,std_curve.iloc[0])]
    std_curve = std_curve.iloc[1:].reset_index(drop=True)

    # melt columns with sparse data into one column
    meltcols = ['Blank (53) MFI', 'HbOHA (27) MFI', 'Pertussis (43) MFI', 'HepB (18) MFI',
           'Tetanus (65) MFI', 'Diptheria (25) MFI', 'RSV-F, A2 (29) MFI',
           'Rubella (52) MFI']
    non_meltcols = list(set(std_curve.columns).difference(meltcols))
    std_curve = std_curve.melt(id_vars = non_meltcols, value_vars = meltcols, var_name='antigen_string', value_name='result_mfi')
    std_curve = std_curve.loc[std_curve.result_mfi.notna()]

    # parse out antigen vs beadset num
    antigen_string = std_curve.antigen_string.str.split('(', expand=True)
    antigen_string.columns = ['antigen', 'beadset_num']
    antigen_string.antigen = antigen_string.antigen.str.strip()
    antigen_string.beadset_num = antigen_string.beadset_num.str.split(")", expand=True)[0].astype(int)
    std_curve = pd.concat([std_curve, antigen_string], axis=1)

    # pull out columns of interest
    std_curve = std_curve.rename(columns={
        'Expected Concentration x Sample Dilution (50)':'expected_concentration',
        'Description':'description_standard'
    })
    std_curve = std_curve[['description_standard', 'antigen','beadset_num','result_mfi','expected_concentration']]

    # add metadata
    std_curve['specrole'] = 'Standard'
    std_curve['MFI Upload Date'] = '2025-03-21'
    std_curve['Concentration Upload Date'] = '2025-03-21'
    std_curve = std_curve.merge(antigen_to_concentration_units, on='antigen', how='left')

    std_curve.loc[std_curve.antigen=='Blank', 'expected_concentration'] = np.nan

    data['specrole'] = 'Sample'
    data = pd.concat([data, std_curve])

    # loq data hard coded from green rows of Antigen Standards tab
    LLOQ = {
        'Diptheria': 0.013717421,
        'HbOHA': 0.045724737,
        'HepB': np.nan,
        'Pertussis': 0.685871056,
        'RSV-F, A2': 0.076207895,
        'Rubella': 1.143118427,
        'Tetanus': 0.022862369,
    }

    ULOQ = {
        'Diptheria': 3.333333333,
        'HbOHA': 33.33333333,
        'HepB': np.nan,
        'Pertussis': 55.55555556,
        'RSV-F, A2': 6.172839506,
        'Rubella': 277.7777778,
        'Tetanus': 1.851851852,
    }

    data["LLOQ_concentration"] = data.antigen.map(LLOQ)
    data["ULOQ_concentration"] = data.antigen.map(ULOQ)

    data.columns = [i.lower().replace(" ","_") for i in data.columns]

    # merge on ldms / standard processing ------------------------------------##
    ldms_feeddir = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN135/specimens/ldms_feed/'
    ldms_path = ldms_feeddir + np.sort(os.listdir(ldms_feeddir))[-1]

    ldms = pd.read_csv(
        ldms_path,
        usecols=constants.STANDARD_COLS,
        dtype=constants.LDMS_DTYPE_MAP
    )

    ldms = ldms.loc[ldms.guspec.isin(data.guspec)]

    md = {
        'network': 'HVTN',
        'upload_lab_id': 'P5',
        'assay_lab_name': 'Fouda Lab (Cornell)',
        'instrument': 'Bio-Rad BioPlex 200',
        'lab_software_version': 'BioPlex 6.2',
        'assay_type': 'BAMA',
        'assay_subtype': 'PVMA',
        'specrole': 'Sample',
        'result_mfi_units': 'RLU'
    }

    outputs = sdmc.standard_processing(
        input_data=data,
        input_data_path=input_data_path,
        guspec_col="guspec",
        network="HVTN",
        metadata_dict=md,
        ldms=ldms
    )

    # this is now a bug in the sdmc-tools adhoc processing pacakge. need to fix
    outputs.loc[outputs.guspec.isna(),'specrole'] = 'Standard'

    drop_cols = ['assay_lab', 'der']
    outputs = outputs.drop(columns=drop_cols)
    outputs = outputs.rename(columns={'mfi_<_100':'mfi<100'})

    # rename to match the ELISA'd HepB data
    outputs = outputs.rename(columns={'sample_dilution':'dilution'})

    # excel turned this into a date; fix.
    outputs.dilution = outputs.dilution.astype(str).str[1:5]

    reorder = [
        'network',
        'protocol',
        'specrole',
        'description_standard',
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
        'instrument',
        'assay_date',
        'assay_tech',
        'internal_qc',
        'lab_internal_sample_id_20240822_upload',
        'lab_internal_sample_id_20250321_upload',
        'lab_software_version',
        'mfi<100',
        'plate_number',
        'antigen',
        'beadset_num',
        'dilution',
        'result_mfi',
        'result_mfi_units',
        'result_concentration',
        'expected_concentration',
        'result_concentration_units',
        'lloq_concentration',
        'uloq_concentration',
        'mfi_upload_date',
        'concentration_upload_date',
        'sdmc_processing_datetime',
        'sdmc_data_receipt_datetime',
        'input_file_name'
    ]
    set(outputs.columns).symmetric_difference(reorder)
    outputs = outputs[reorder]

    # save to .txt file ------------------------------------------------------##
    savedir = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN135/assays/bama_pvma/misc_files/data_processing/'
    today = datetime.date.today().isoformat()

    outputs.to_csv(savedir + f"HVTN135_PVMA_FOUDA_processed_{today}.txt", index=False, sep="\t")

    # save pivot summary -----------------------------------------------------##
    pivot_summary = outputs.copy()
    pivot_summary = pivot_summary.replace("See ELISA Data",np.nan)
    pivot_summary[['ptid','visitno']] = pivot_summary[['ptid','visitno']].fillna("NA")
    pivot_summary = pd.pivot_table(data=pivot_summary,
                   index=['specrole','visitno','ptid'],
                   columns=['antigen'],
                   aggfunc='count',
                   fill_value=0,
                  )[['result_mfi','result_concentration']]

    pivot_summary.to_excel(savedir + "HVTN135_PVMA_FOUDA_pivot_summary.xlsx")

def checks():
    savedir = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN135/assays/bama_pvma/misc_files/data_processing/'
    files = os.listdir(savedir)
    fname = [i for i in files if 'HVTN135_PVMA_FOUDA_processed' in i][0]
    output = pd.read_csv(savedir + fname, sep="\t")

    manifest_path = "/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN135/assays/bama_pvma/misc_files/HVTN 135 v12 Shipment to Cornell_Fouda.xlsx"
    manifest = pd.read_excel(manifest_path, sheet_name="Sheet2")
    manifest_mismatch = set(outputs.guspec.dropna()).symmetric_difference(manifest.GLOBAL_ID)
    if len(manifest_mismatch) > 0:
        raise Exception("Manifest doesn't match data received")

    original_data_path = "/trials/vaccine/p135/s001/qdata/LabData/PVMA_pass-through/preliminary_data/April2024_HVTN 135 Summary Data_V1.xlsx"
    og_data = pd.read_excel(original_data_path)
    og_mismatch = set(outputs.guspec.dropna()).symmetric_difference(og_data.GLOBAL_ID)
    if len(og_mismatch) > 0:
        raise Exception("Samples don't match prelimiary data")

if __name__=="__main__":
    main()
    checks()
