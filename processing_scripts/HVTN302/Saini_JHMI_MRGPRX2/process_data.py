import sdmc_tools.constants as constants
import sdmc_tools.process as sdmc

import pandas as pd
import numpy as np

import os
import datetime

pd.options.mode.chained_assignment = None  # default='warn'

def main():
    # read in data
    datadir = '/trials/vaccine/p302/s001/qdata/LabData/AE_assays_pass-through/MGPRX2/'
    fnames = [i for i in os.listdir(datadir) if i[:4]=='HVTN']

    def get_ppt_list(fname):
        """
        given one of the input data filenames,
        return a list of the RUNS contained in the data as ints
        """
        ppts = []

        idx, l = 0, 0
        remainder = fname[idx+l:]
        while "#" in remainder:
            idx = remainder.index("#")
            l = remainder[idx:].index("-")
            ppts += [int(remainder[idx+1:idx+l])]
            remainder = remainder[idx + l:]
        return ppts

    def get_idx_of_Average(df):
        """
        Find location (i,j) in the sheet where the word "Average"
        located (this is where romi put the final data)
        """
        np_version = np.array(df)
        if 'average' in np_version:
            return np.where(np_version=='average')
        elif 'Average' in np_version:
            return np.where(np_version=='Average')
        else:
            return "NONE FOUND"

    def get_subset(f):
        """
        Given input data filename,
        return subset of excel sheet with the 'Average' data table
        """
        df = pd.read_excel(datadir + f, engine='calamine')
        i, j = get_idx_of_Average(df)
        return df.iloc[i[0]:,j[0] - 1: j[0] + 2]

    all_data = {}
    for f in fnames:
        all_data[f] = get_subset(f)

    # pull table of number of visits per "RUN" id
    sample_mapping_file_path = '/trials/vaccine/p302/s001/qdata/LabData/AE_assays_pass-through/MGPRX2/Summary IHVTN samples MRGPRX2.xlsx'
    sample_mapping = pd.read_excel(sample_mapping_file_path)
    sample_mapping = sample_mapping.iloc[3:].dropna(how='all')
    sample_mapping.columns = sample_mapping.iloc[0]
    sample_mapping = sample_mapping.iloc[1:]
    nruns = sample_mapping[['# visits', 'RUNS']]

    # subset to only desired data in sheet
    def grab_desired_data(f):
        df = all_data[f]
        # set WT vs KO as column names
        df.columns = df.iloc[1]

        # drop header columns
        df = df.iloc[2:]

        # keep only sample data + 2 rows controls
        ppts = get_ppt_list(f)
        N_samples = nruns.loc[nruns.RUNS.isin(ppts),'# visits'].unique()[0]
        k = len(ppts)*N_samples + 2
        df = df.iloc[:k]

        return df

    pdata = {}
    for f in fnames:
        pdata[f] = grab_desired_data(f)

    # turn the positive + maximal control rows into columns
    def flip_control_columns_to_wide(f):
        rename_map = {
        'KO_ 5ug/ml C48/80': 'KO_C48/80',
        'KO_1% Tween20': 'KO_1pct_Tween20',
        'KO_C48/80': 'KO_C48/80',
        'KO_C48/80 ': 'KO_C48/80',
        'KO_Tween 20': 'KO_1pct_Tween20',
        'KO_Tween20  1%': 'KO_1pct_Tween20',
        'KO_Tween20 1%': 'KO_1pct_Tween20',
        'WT_ 5ug/ml C48/80': 'WT_C48/80',
        'WT_1% Tween20': 'WT_1pct_Tween20',
        'WT_C48/80': 'WT_C48/80',
        'WT_C48/80 ': 'WT_C48/80',
        'WT_Tween 20': 'WT_1pct_Tween20',
        'WT_Tween20  1%': 'WT_1pct_Tween20',
        'WT_Tween20 1%': 'WT_1pct_Tween20'
        }

        df = pdata[f]
        df = df.rename(columns={ df.columns[0]: "jhmi_sample_id" })

        s = df.iloc[-2:]
        s = s.melt(id_vars=['jhmi_sample_id'], value_vars=['WT','KO'], var_name='wt_vs_ko', value_name='result')
        s['sample_type'] = s.wt_vs_ko + "_" + s.jhmi_sample_id
        s['sample_type'] = s['sample_type'].map(rename_map)
        s = s.set_index('sample_type')[['result']].T
        t = df.iloc[:-2]
        return t.merge(s, how='cross')

    pdata2 = {}
    for f in fnames:
        pdata2[f] = flip_control_columns_to_wide(f)

    df = pd.concat([pdata2[f] for f in fnames], axis = 0)

    # merge on guspec col
    hvtn_sample_mapping_path = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN302/assays/AE_assays/MRGPRX2/misc_files/confirm_sample_mapping_2024_12_02.xlsx'
    hvtn_sample_mapping = pd.read_excel(hvtn_sample_mapping_path, skiprows=1)

    hvtn_sample_mapping = hvtn_sample_mapping.rename(columns={'HVTN_sample_id':'guspec',
                                                              'JHMI_sample_id':'jhmi_sample_id'
                                                             })
    df = df.merge(hvtn_sample_mapping[['guspec','jhmi_sample_id']], on=['jhmi_sample_id'], how='outer')

    ## STANDARD PROCESSING ---------------------------------------------------##
    ldms = pd.read_csv(constants.LDMS_PATH_HVTN, usecols=constants.STANDARD_COLS, dtype=constants.LDMS_DTYPE_MAP)
    ldms = ldms.loc[ldms.lstudy==302.]
    ldms = ldms.loc[ldms.guspec.isin(df.guspec)]

    df = df.rename(columns={
        'WT': 'result_WT',
        'KO': 'result_KO',
        'WT_C48/80': 'result_WT_C48/80',
        'WT_1pct_Tween20': 'result_WT_1pct_Tween20',
        'KO_C48/80': 'result_KO_C48/80',
        'KO_1pct_Tween20': 'result_KO_1pct_Tween20'
    })

    input_filepath_example = '/trials/vaccine/p302/s001/qdata/LabData/AE_assays_pass-through/MGPRX2/HVTN302 MRGPRX2 Run-#1-#2-S.xlsx'
    md = {
        'network': 'HVTN',
        'specrole': 'Sample',
        'upload_lab_id': 'X3',
        'assay_lab_name': 'Saini Lab (JHMI)',
        'assay_type': 'MRGPRX2',
        'instrument': 'Molecular Devices FlexStation 3',
        'lab_software_version': 'SoftMax Pro 7.02',
        'result_units': 'Percent'
    }

    outputs = sdmc.standard_processing(
                input_data=df,
                input_data_path=input_filepath_example,
                additional_input_paths={'jhmi_sample_mapping_file':sample_mapping_file_path,
                                        'hvtn_sample_mapping_file':hvtn_sample_mapping_path},
                guspec_col="guspec",
                network="HVTN",
                metadata_dict=md,
                ldms=ldms
    )

    # correct input_file_name column with ppt-specific input file
    ppts_to_files = {}
    for f in fnames:
        ppts = get_ppt_list(f)
        for p in ppts:
            ppts_to_files[p] = f
    outputs['input_file_name'] = outputs.jhmi_sample_id.str.replace("-","_").str.split("_", expand=True)[0].astype(int).map(ppts_to_files)

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
        'result_ko',
        'result_ko_1pct_tween20',
        'result_ko_c48/80',
        'result_wt',
        'result_wt_1pct_tween20',
        'result_wt_c48/80',
        'result_units',
        'jhmi_sample_id',
        'jhmi_sample_mapping_file',
        'hvtn_sample_mapping_file',
        'sdmc_processing_datetime',
        'sdmc_data_receipt_datetime',
        'input_file_name'
    ]

    mismatch = set(outputs.columns).symmetric_difference(reorder)
    if len(mismatch) > 0:
        raise Warning(f"Trying to add or remove columns! Problem cols: {mismatch}")
    outputs = outputs[reorder]

    ## save to .txt
    savedir = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN302/assays/AE_assays/MRGPRX2/misc_files/data_processing/'
    today = datetime.date.today().isoformat()
    fname = f'HVTN302_Saini_JHMI_MRGPRX2_Processed_{today}.txt'

    outputs.to_csv(savedir + fname, index=False, sep="\t")

    ## save pivot summary of ppts/visits in dataset
    summary = pd.pivot_table(outputs, index=['ptid'], columns='visitno', aggfunc='count', fill_value=0)[['result_wt']].droplevel(level=0, axis=1)
    summary = summary[['2.0', '3.0', '4.0', '6.0', '7.0', '8.0', '10.0', '12.0', '15.0',]]
    summary.to_excel(savedir + "HVTN302_Saini_JHMI_MRGPRX2_summary.xlsx")

    # check for completeness against shipping manifest
    shipping_manifest_path = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN302/assays/AE_assays/MRGPRX2/misc_files/johns_hopkins_mrgprx_shipping_manifest.txt'
    shipping_manifest = pd.read_csv(shipping_manifest_path, sep="\t")
    mismatch = set(outputs.guspec).symmetric_difference(shipping_manifest.GLOBAL_ID)
    if len(mismatch) > 0:
        raise Warning(f"data doesn't match shipping manifest!")

if __name__=="__main__":
    main()
