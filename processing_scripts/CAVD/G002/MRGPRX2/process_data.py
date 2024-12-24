## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 2024-12-23
# Purpose: Process G002 MRGPRX2 data from Saini Lab
## ---------------------------------------------------------------------------##
import pandas as pd
import numpy as np

import os
import datetime

def main():
    ## pull in data and merge together -------------------------------------- ##
    datadir = '/networks/cavd/VDCs/Schief/Schief_856-G002/SkinReactions/data/MRGPX2_Saini_JHU/'
    fnames = [i for i in os.listdir(datadir) if i[:4]=='2024']

    def get_idx_of_Average(df):
        """
        Find location (i,j) in the sheet where the word "Average"
        located (this is where romi put the final data)
        """
        np_version = np.array(df)
        if 'ave' in np_version:
            return np.where(np_version=='ave')
        elif 'Ave' in np_version:
            return np.where(np_version=='Ave')
        else:
            return "NONE FOUND"

    def get_subset(f):
        d = pd.read_excel(f, engine='calamine')

        i, j = get_idx_of_Average(d)
        i, j = i[0], j[0]

        return d.iloc[i:,j-2:j+2].dropna(how='all')

    all_data = {}
    for f in fnames:
        all_data[f] = get_subset(datadir + f)

    def grab_desired_data(f):
        d = all_data[f]
        d = d.reset_index(drop=True)
        d.columns = d.iloc[1]
        d = d.iloc[2:]
        d.columns = ['jhmi_sample_id_a', 'jhmi_sample_id_b', 'wt', 'ko']

        s = d.iloc[-2:]
        d = d.iloc[:-2]

        s = s.melt(id_vars=['jhmi_sample_id_b'], value_vars=['wt','ko'])
        s.columns = ['control', 'wt_vs_ko', 'result']
        s['sample_type'] = s.wt_vs_ko + "_" + s.control.str.rstrip().str.replace("%","pct").replace(" ","_")

        s = s.set_index('sample_type')[['result']].T
        d = d.merge(s, how='cross')
        d['input_file_name'] = f
        return d

    ## merge on metadata from sample sheet ---------------------------------- ##
    df = pd.concat([grab_desired_data(f) for f in fnames])
    df[['ptid','visitno']] = df.jhmi_sample_id_a.str.split("_", expand=True)

    sample_sheet_fname = 'Summary IAVI samples 12.2024.xlsx'
    sample_sheet = pd.read_excel(datadir + sample_sheet_fname)
    sample_sheet.columns = sample_sheet.iloc[3]
    sample_sheet = sample_sheet.iloc[4:].dropna(how='all')
    sample_sheet = sample_sheet.reset_index(drop=True)

    for i in range(len(sample_sheet)):
        if pd.isnull(sample_sheet.Site.iloc[i]):
            sample_sheet.Site.iloc[i] = sample_sheet.Site.iloc[i-1]

    box_map = sample_sheet.dropna().set_index('Site')['SAMPLE'].to_dict()
    sample_sheet.SAMPLE = sample_sheet.Site.map(box_map)
    sample_sheet.columns = ['Site', 'ptid', 'visit_count','volume (mL)', 'Box', 'completed']

    df['potential_error'] = False
    df['error_detail'] = 'N/A'

    # FIXING TYPOS -- EXPECTING TO UPDATE AFTER HEARING FROM ROMI ----------- ##
    df.loc[df.ptid.isin(['G259002', 'G259005']),'potential_error'] = True
    df.loc[df.ptid.isin(['G259002']),'error_detail'] = "Lab submitted ptid 'G259002'"
    df.loc[df.ptid.isin(['G259005']),'error_detail'] = "Lab submitted ptid 'G259005'"

    df.loc[df.ptid.isin(['G259002', 'G259005']),'ptid'] = df.loc[df.ptid.isin(['G259002', 'G259005'])].ptid.map({'G259002':'G00259002', 'G259005':'G00259005'})
    df = df.merge(sample_sheet[['Site','Box','ptid']], on='ptid', how='left')

    # # make sure sample sheet boxes map individual sheet boxes
    # df.loc[df.jhmi_sample_id_b.str[0]!=df.Box.str[-1]]

    # # make sure count of datapoints per ptid matches expectations from sample sheet
    # t = df.groupby('ptid').count()[['jhmi_sample_id_a']].reset_index()
    # compare = sample_sheet.merge(t, on='ptid', how='outer')
    # compare.loc[compare['visit_count'].astype(float)!=compare.jhmi_sample_id_a.astype(float)]

    ## merge on metadata from primary sample sheet -------------------------- ##
    sample_sheet_path = '/networks/cavd/VDCs/Schief/Schief_856-G002/SkinReactions/data/MRGPX2_Saini_JHU/Summary IAVI samples 04.2024.xlsx'
    sample_ref = pd.read_excel(sample_sheet_path, sheet_name="All samples listed", skiprows=114)

    # grab each set (box) from excel sheet
    sample_ref1 = sample_ref.iloc[:29]
    sample_ref2 = sample_ref.iloc[29:95].reset_index(drop=True)
    sample_ref3 = sample_ref.iloc[95:120].reset_index(drop=True)
    sample_ref4 = sample_ref.iloc[120:].reset_index(drop=True)

    sample_ref2.columns = sample_ref2.iloc[0]
    sample_ref2 = sample_ref2.iloc[1:]

    sample_ref3.columns = sample_ref3.iloc[0]
    sample_ref3 = sample_ref3.iloc[1:]

    sample_ref4.columns = sample_ref4.iloc[0]
    sample_ref4 = sample_ref4.iloc[1:]

    # rename columns to standardize
    sample_ref1 = sample_ref1.rename(columns={'LIMS ID':'sample_id',
                                              'PTID':'ptid',
                                              'VISIT#':'visitno',
                                              'Collection date/time':'drawdt'
                                             })

    sample_ref2 = sample_ref2.rename(columns={
        'Unique Vial ID':'sample_id',
        'PTID':'ptid',
        'Visit':'visitno',
        'Date Collected': 'drawdt'
    })

    sample_ref3 = sample_ref3.rename(columns={
        'Global Spec ID':'sample_id',
        'ID1':'ptid',
        'Visit':'visitno',
        'Specimen Date': 'drawdt'
    })

    sample_ref4 = sample_ref4.rename(columns={
        'Unique Aliquot ID':'sample_id',
        'PTID':'ptid',
        'Visit':'visitno',
        'Date Collected':'drawdt'
    })

    #merge together
    usecols = ['sample_id','ptid','visitno','drawdt']
    ref_sheet = pd.concat([sample_ref1[usecols], sample_ref2[usecols], sample_ref3[usecols], sample_ref4[usecols]])

    # checking for agreement between sheets
    df.visitno = df.visitno.astype(int)
    ref_sheet.visitno = ref_sheet.visitno.astype(int)

    ref_sheet.ptid = ref_sheet.ptid.str.replace("-","").str.strip()
    df.ptid = df.ptid.astype(str).str.strip()

    tmp = df[['ptid','visitno']].merge(ref_sheet, on=['ptid','visitno'], how='outer')

    # CORRECT TYPOS -- EXPECTING TO UPDATE AFTER HEARING FROM ROMI ---------- ##
    df.loc[(df.ptid=="G00258016") & (df.visitno==160), 'potential_error'] = True
    df.loc[(df.ptid=="G00258016") & (df.visitno==160), 'error_detail'] = "Lab submitted ptids 160 and 161 for this row. 161 appears to be correct; this was an ad hoc visit likely due to insufficient blood drawn during visit 160 on 2022-07-13"
    df.loc[(df.ptid=="G00258016") & (df.visitno==160), 'visitno'] = 161

    df.loc[(df.ptid=="G00258025") & (df.visitno==550), 'potential_error'] = True
    df.loc[(df.ptid=="G00258025") & (df.visitno==550), 'error_detail'] = "Lab submitted visitnos 550 and 540 for this row"
    df.loc[(df.ptid=="G00258025") & (df.visitno==550), 'visitno'] = 540

    df.loc[(df.ptid=="G00258026") & (df.visitno==550), 'potential_error'] = True
    df.loc[(df.ptid=="G00258026") & (df.visitno==550), 'error_detail'] = "Lab submitted visitnos 550 and 540 for this row"
    df.loc[(df.ptid=="G00258026") & (df.visitno==550), 'visitno'] = 540

    df = df.merge(ref_sheet, on=['ptid','visitno'], how = 'left')

    ## final formatting + metadata ------------------------------------------ ##
    df.columns = [i.lower().replace(" ","_") for i in df.columns]

    receipt_time = datetime.datetime.fromtimestamp(os.path.getmtime(datadir + sample_sheet_fname)).replace(microsecond=0).isoformat()
    sdmc_processing_timestamp = datetime.datetime.now().replace(microsecond=0).isoformat()

    metadata = {
        'network': ['CAVD'],
        'protocol': ['G002'],
        'specrole': ['Sample'],
        'upload_lab_id': ['X3'],
        'assay_lab_name': ['Saini Lab (JHMI)'],
        'assay_type': ['MRGPRX2'],
        'instrument': ['Molecular Devices FlexStation 3'],
        'lab_software_version': ['SoftMax Pro 7.02'],
        'lab_result_units': ['Percent'],
        'jhmi_sample_mapping_file': ['Summary IAVI samples 04.2024.xlsx'],
        'jhmi_secondary_mapping_file': [sample_sheet_fname],
        'sdmc_data_receipt_datetime': [receipt_time],
        'sdmc_processing_datetime': [sdmc_processing_timestamp],
    }
    df = df.merge(pd.DataFrame(metadata), how='cross')

    rename = {
        'wt': 'lab_result_wt',
        'ko': 'lab_result_ko',
        'wt_1pct_tween_20': 'lab_result_wt_1pct_tween20',
        'ko_1pct_tween_20': 'lab_result_ko_1pct_tween20',
        'wt_5ug/ml_c48/80':'pc_wt_c48/80',
        'ko_5ug/ml_c48/80':'pc_ko_c48/80',
    }
    df = df.rename(columns=rename)

    reorder = [
        'network',
        'protocol',
        'specrole',
        'sample_id',
        'ptid',
        'visitno',
        'drawdt',
        'upload_lab_id',
        'assay_lab_name',
        'assay_type',
        'instrument',
        'lab_software_version',
        'site',
        'box',
        'lab_result_ko',
        'lab_result_ko_1pct_tween20',
        'lab_result_wt',
        'lab_result_wt_1pct_tween20',
        'lab_result_units',
        'pc_ko_c48/80',
        'pc_wt_c48/80',
        'jhmi_sample_id_a',
        'jhmi_sample_id_b',
        'jhmi_sample_mapping_file',
        'jhmi_secondary_mapping_file',
        'sdmc_data_receipt_datetime',
        'sdmc_processing_datetime',
        'input_file_name',
        'potential_error',
        'error_detail',
    ]

    mismatch = set(df.columns).symmetric_difference(reorder)
    if len(mismatch) > 0:
        raise Warning(f"Trying to add or remove columns: {mismatch}")

    df = df[reorder]

    ## save ----------------------------------------------------------------- ##
    savedir = '/networks/vtn/lab/SDMC_labscience/studies/VISC/G002/MRGPRX2/data_processing/'
    today = datetime.date.today().isoformat()
    fname = f'PROVISIONAL_CAVD_G002_MRGPRX2_Processed_{today}.txt'

    df.to_csv(savedir + fname, sep="\t", index=False)

    df_calculated = df.copy()

    df_calculated['result_wt_normalized'] = df_calculated.lab_result_wt/df_calculated.lab_result_wt_1pct_tween20
    df_calculated['result_ko_normalized'] = df_calculated.lab_result_ko/df_calculated.lab_result_ko_1pct_tween20
    df_calculated['result_delta_pct_maximal_release'] = df_calculated.result_wt_normalized - df_calculated.result_ko_normalized
    df_calculated['result_units'] = "Percent (Normalized)"

    reorder_calculated = [
        'network',
        'protocol',
        'specrole',
        'sample_id',
        'ptid',
        'visitno',
        'drawdt',
        'upload_lab_id',
        'assay_lab_name',
        'assay_type',
        'instrument',
        'lab_software_version',
        'site',
        'box',
        'lab_result_ko',
        'lab_result_ko_1pct_tween20',
        'lab_result_wt',
        'lab_result_wt_1pct_tween20',
        'lab_result_units',
        'result_ko_normalized',
        'result_wt_normalized',
        'result_delta_pct_maximal_release',
        'result_units',
        'pc_ko_c48/80',
        'pc_wt_c48/80',
        'jhmi_sample_id_a',
        'jhmi_sample_id_b',
        'jhmi_sample_mapping_file',
        'jhmi_secondary_mapping_file',
        'sdmc_data_receipt_datetime',
        'sdmc_processing_datetime',
        'input_file_name',
        'potential_error',
        'error_detail',
    ]

    mismatch = set(df_calculated.columns).symmetric_difference(reorder_calculated)
    if len(mismatch) > 0:
        raise Warning(f"Trying to add or remove columns: {mismatch}")

    df_calculated = df_calculated[reorder_calculated]

    today = datetime.date.today().isoformat()
    fname_calculated = f'PROVISIONAL_CAVD_G002_MRGPRX2_calculated_delta_pct_release_{today}.txt'

    df_calculated.to_csv(savedir + fname_calculated, sep="\t", index=False)

    ## pivot summary -------------------------------------------------------- ##
    # pivot_summary = pd.pivot_table(df, index='ptid', columns='visitno', aggfunc='count', fill_value=0)['wt']
    # pivot_summary.to_excel(savedir + 'CAVD_G002_MRGPRX_ptid_visit_summary.xlsx')

    pivot_summary = pd.pivot_table(df, index='ptid', columns='visitno', aggfunc='count', fill_value=0)['lab_result_wt']
    pivot_summary.to_excel(savedir + 'CAVD_G002_MRGPRX_ptid_visit_summary_2024_12_23.xlsx')

if __name__=="__main__":
    main()
