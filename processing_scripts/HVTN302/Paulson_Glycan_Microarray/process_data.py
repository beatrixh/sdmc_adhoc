## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 01.21.2025
# Purpose:
## ---------------------------------------------------------------------------##
import pandas as pd
import numpy as np
import os
import datetime
import sdmc_tools.process as sdmc
import sdmc_tools.constants as constants

## ---------------------------------------------------------------------------##
def get_data(path, channel):
    nrows = 35
    data = pd.read_csv(path, encoding='iso-8859-1', sep="\t", skiprows=nrows)
    while 'Block' not in data.columns and nrows > 0:
        nrows -= 1
        data = pd.read_csv(path, encoding='iso-8859-1', sep="\t", skiprows=nrows)
    data['ptid'] = path.split("/")[-1].split("-")[-1].split("_")[0]
    data.ptid = data.ptid.astype(int)
    data['filter_pass'] = path.split("/")[-1].split(".")[0].split("_")[-1].strip()
    data['filename'] = path.split("/")[-1]
    data['signal_channel'] = channel
    data = data.rename(columns={'Name':'spot_name',
                                f'F{channel} Mean': f'mean_signal',
                                f'B{channel} Mean': f'mean_background_signal',
                                'Block':'block'})
    usecols = ['ptid', 'spot_name', 'block', 'filter_pass', 'signal_channel', 'mean_signal', 'mean_background_signal','filename']
    data = data.loc[~data.spot_name.isin(["empty","Empty","DyeSpot"]),usecols]

    return data


def pull_and_merge_inputs(datadir, channel, sample_mapping_path, glycan_info_path, filter_pass):
    data = pd.DataFrame()
    files = os.listdir(datadir)
    files = [f for f in files if f[-3:]=="txt"]
    for file in files:
        data = pd.concat([data, get_data(datadir + file, channel)])

    # from the lab: we only want the "Hi" data for IgG/IgE, "Lo" for IgM
    data = data.loc[data.filter_pass==filter_pass]

    metadata = pd.read_excel(sample_mapping_path)
    metadata.columns = [i.lower().replace(" ","_") for i in metadata.columns]
    metadata = metadata.rename(columns = {'original_id':'guspec', 'subject_id':'ptid', 'box_':'box'})
    metadata['visitno'] = metadata['patient-visit'].str.split("-V", expand=True)[1].astype(int)

    def label_timepoint(ptid, visitno):
        visits = metadata.loc[(metadata.ptid==ptid)].visitno.tolist()
        return visits.index(visitno)

    metadata['timepoint'] = metadata.apply(lambda x: label_timepoint(x['ptid'], x['visitno']), axis=1)
    block_to_timepoint = {**{i:0 for i in range(13)}, **{i:1 for i in range(13,25)}, **{i:2 for i in range(25,37)}}
    data["timepoint"] = data.block.map(block_to_timepoint)
    data = data.merge(metadata[['ptid','timepoint','guspec']], on=['ptid','timepoint'], how='left')

    # dropping merge column
    data = data.drop(columns='timepoint')

    # dropping rows for which don't have data
    data = data.loc[data.guspec.notna()]

    # ## read glycan info ------------------------------------------------------##
    glyc_info = pd.read_excel(glycan_info_path)
    glyc_info = glyc_info.rename(columns={'Sample':'s', 'M-Number':'glycan_m_number', 'Structure':'glycan_structure'})

    def replace_chars(s):
        return s.replace("α","[alpha]").replace("β","[beta]").replace("\xa0"," ").replace('Α','A')

    glyc_info.loc[glyc_info.glycan_structure.isna(),'glycan_structure'] = 'Not provided'
    glyc_info.glycan_structure = glyc_info.glycan_structure.apply(lambda x: replace_chars(x))

    # merge on M-Number from glyc_info
    tmp = data[['spot_name']].drop_duplicates()
    tmp['s'] = tmp.spot_name.str.partition("-", expand=True)[0]
    tmp.s = tmp.s.astype(int)
    tmp['data_glycan'] = tmp.spot_name.str.partition("-", expand=True)[2]
    tmp = glyc_info.merge(tmp, on='s', how='outer')

    Name_to_glycan_m_number_map = {data:ref for (data,ref) in zip(tmp.spot_name, tmp.glycan_m_number)}
    data['glycan_m_number'] = data.spot_name.map(Name_to_glycan_m_number_map)

    data = data.merge(glyc_info[['glycan_m_number','glycan_structure']], on='glycan_m_number', how='left')

    ## final pre-processing --------------------------------------------------##

    # add primary result column
    # there's a typo here but maintaining it because stats has already coded off it
    data[f'background_subtraced_mean_signal'] = data.mean_signal - data.mean_background_signal
    data['background_subtraced_mean_signal'] = data['background_subtraced_mean_signal'].round(2)

    return data

def main():
    ## read in ldms ----------------------------------------------------------##
    ldms = pd.read_csv(constants.LDMS_PATH_HVTN, usecols=constants.STANDARD_COLS, dtype=constants.LDMS_DTYPE_MAP)
    ldms = ldms.loc[ldms.lstudy==302.]

    ## hand-entered data -----------------------------------------------------##
    metadata_dict = {
        "network": "HVTN",
        "upload_lab_id": "P4",
        "assay_lab_name": "Paulson (Scripps)",
        "assay_type": "Glycan Microarray",
        "assay_precision": "Semi-Quantitative",
        "assay_details":"Custom glycan layout",
        "instrument": "Innoscan 1100AL",
        "lab_software_version": "Mapix",
        "specrole": "Sample",
        "glycan_mapping_file": "Scheif-GlycanArrayList-Final.xlsx",
        "sample_mapping_file": "HVTN302-PatientSample-GlycanAnalysis.xlsx"
    }

    ## IgE -------------------------------------------------------------------##
    ige_datadir = '/trials/vaccine/p302/s001/qdata/LabData/AE_assays_pass-through/Glycan_array/IgE/lab_processed_data/'
    sample_mapping_path = "/trials/vaccine/p302/s001/qdata/LabData/AE_assays_pass-through/Glycan_array/IgE/lab_processed_data/HVTN302-PatientSample-GlycanAnalysis.xlsx"
    glycan_info_path = "/trials/vaccine/p302/s001/qdata/LabData/AE_assays_pass-through/Glycan_array/metadata_approved_by_lab/Scheif-GlycanArrayList-Final.xlsx"

    ige = pull_and_merge_inputs(ige_datadir, 532, sample_mapping_path, glycan_info_path, "Hi")
    ige['isotype'] = 'IgE'

    files = os.listdir(ige_datadir)
    files = [f for f in files if f[-3:]=="txt"]

    ige_outputs = sdmc.standard_processing(
        input_data=ige.drop(columns='ptid'),
        input_data_path=ige_datadir + files[0],
        guspec_col="guspec",
        network="HVTN",
        metadata_dict=metadata_dict,
        ldms=ldms.loc[ldms.guspec.isin(ige.guspec)]
    )

    ## IgG -------------------------------------------------------------------##
    igg_datadir = '/trials/vaccine/p302/s001/qdata/LabData/AE_assays_pass-through/Glycan_array/IgG/lab_processed_data/'
    sample_mapping_path = "/trials/vaccine/p302/s001/qdata/LabData/AE_assays_pass-through/Glycan_array/IgG/lab_processed_data/metadata/HVTN302-PatientSample-GlycanAnalysis.xlsx"
    glycan_info_path = "/trials/vaccine/p302/s001/qdata/LabData/AE_assays_pass-through/Glycan_array/metadata_approved_by_lab/Scheif-GlycanArrayList-Final.xlsx"

    igg = pull_and_merge_inputs(igg_datadir, 488, sample_mapping_path, glycan_info_path, "Hi")
    igg['isotype'] = 'IgG'

    files = os.listdir(igg_datadir)
    files = [f for f in files if f[-3:]=="txt"]

    igg_outputs = sdmc.standard_processing(
        input_data=igg.drop(columns='ptid'),
        input_data_path=igg_datadir + files[0],
        guspec_col="guspec",
        network="HVTN",
        metadata_dict=metadata_dict,
        ldms=ldms.loc[ldms.guspec.isin(igg.guspec)]
    )

    ## IgM old -------------------------------------------------------------------##
    igm_datadir = '/trials/vaccine/p302/s001/qdata/LabData/AE_assays_pass-through/Glycan_array/IgM/lab_processed_data/10_01_2024/'
    sample_mapping_path = '/trials/vaccine/p302/s001/qdata/LabData/AE_assays_pass-through/Glycan_array/IgM/HVTN302-PatientSample-GlycanAnalysis.xlsx'
    glycan_info_path = "/trials/vaccine/p302/s001/qdata/LabData/AE_assays_pass-through/Glycan_array/metadata_approved_by_lab/Scheif-GlycanArrayList-Final.xlsx"

    igm = pull_and_merge_inputs(igm_datadir, 532, sample_mapping_path, glycan_info_path, "Lo")
    igm['isotype'] = 'IgM'

    files = os.listdir(igm_datadir)
    files = [f for f in files if f[-3:]=="txt"]

    igm_outputs = sdmc.standard_processing(
            input_data=igm.drop(columns='ptid'),
            input_data_path=igm_datadir + files[0],
            guspec_col="guspec",
            network="HVTN",
            metadata_dict=metadata_dict,
            ldms=ldms.loc[ldms.guspec.isin(igm.guspec)]
    )

    ## IgM new -------------------------------------------------------------------##
    igm_datadir_new = '/trials/vaccine/p302/s001/qdata/LabData/AE_assays_pass-through/Glycan_array/IgM/lab_processed_data/11_25_2024/'
    sample_mapping_path = '/trials/vaccine/p302/s001/qdata/LabData/AE_assays_pass-through/Glycan_array/IgM/HVTN302-PatientSample-GlycanAnalysis.xlsx'
    glycan_info_path = "/trials/vaccine/p302/s001/qdata/LabData/AE_assays_pass-through/Glycan_array/metadata_approved_by_lab/Scheif-GlycanArrayList-Final.xlsx"

    igm_new = pull_and_merge_inputs(igm_datadir_new, 488, sample_mapping_path, glycan_info_path, "Lo")
    igm_new['isotype'] = 'IgM'

    files = os.listdir(igm_datadir_new)
    files = [f for f in files if f[-3:]=="txt"]

    igm_new_outputs = sdmc.standard_processing(
            input_data=igm_new.drop(columns='ptid'),
            input_data_path=igm_datadir_new + files[0],
            guspec_col="guspec",
            network="HVTN",
            metadata_dict=metadata_dict,
            ldms=ldms.loc[ldms.guspec.isin(igm_new.guspec)]
    )

    ## combine all -----------------------------------------------------------##

    outputs = pd.concat([ige_outputs, igg_outputs, igm_outputs, igm_new_outputs])

    outputs = outputs.drop(columns="input_file_name")
    outputs = outputs.rename(columns={"filename":"input_file_name"})

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
        'assay_precision',
        'assay_details',
        'spot_name',
        'glycan_m_number',
        'glycan_structure',
        'isotype',
        'signal_channel',
        'mean_signal',
        'mean_background_signal',
        'background_subtraced_mean_signal',
        'filter_pass',
        'block',
        'sample_mapping_file',
        'glycan_mapping_file',
        'input_file_name',
        'sdmc_processing_datetime',
        'sdmc_data_receipt_datetime',
    ]

    mismatch = set(reorder).symmetric_difference(outputs.columns)
    if len(mismatch) > 0:
        raise Exception("Trying to add or remove columns")
    outputs = outputs[reorder]

    ## save ----------------------------------------------------------------------##
    today = datetime.date.today().isoformat()
    fname = f"DRAFT_HVTN302_Glycan_Data_Processed_{today}.txt"
    savedir = "/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN302/assays/AE_assays/glycan_microarray/misc_files/data_processing/"

    outputs.to_csv(savedir + fname, sep="\t", index=False)

    guspec_by_glycan = pd.pivot_table(
                            outputs,
                            index='guspec',
                            columns=['isotype','glycan_m_number'],
                            values='background_subtraced_mean_signal',
                            aggfunc='count',
                            fill_value=0
                            )
    ptid_visit_glycan = pd.pivot_table(
                            outputs,
                            index=['ptid','visitno'],
                            columns=['isotype','glycan_m_number'],
                            aggfunc='count',
                            fill_value=0
                            )
    guspec_by_glycan.to_excel(savedir + "HVTN302_Glycan_Microarray_guspec_summary_count.xlsx")
    ptid_visit_glycan.to_excel(savedir + "HVTN302_Glycan_Microarray_ptid_visitno_summary_count.xlsx")

    shipping_manifest = pd.read_csv('/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN302/assays/AE_assays/glycan_microarray/misc_files/hvtn302_glycan_shipping_manifest_from_nick.txt', sep='\t')
    missings = set(shipping_manifest.GLOBAL_ID).difference(igm_tmp.guspec)
    shipping_manifest.loc[shipping_manifest.GLOBAL_ID.isin(missings),['GLOBAL_ID','PID','VID']].sort_values(by=['PID','VID'])

if __name__=="__main__":
    main()
