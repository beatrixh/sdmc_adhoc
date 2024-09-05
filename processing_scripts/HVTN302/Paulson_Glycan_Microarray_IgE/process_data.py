## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 08/27/2024
# Purpose:  - Concatenate / merge IgE Glycan microarray data
## ---------------------------------------------------------------------------##
import pandas as pd
import numpy as np
import os
import datetime
import sdmc_tools.refactor_process as sdmc
import sdmc_tools.constants as constants
## ---------------------------------------------------------------------------##

def main():
    datadir = '/trials/vaccine/p302/s001/qdata/LabData/AE_assays_pass-through/Glycan_array/IgE/lab_processed_data/'
    files = os.listdir(datadir)
    files = [f for f in files if f[-3:]=="txt"]

    # read in data
    def get_data(path):
        nrows = 35
        data = pd.read_csv(path, encoding='iso-8859-1', sep="\t", skiprows=nrows)
        while 'Block' not in data.columns and nrows > 0:
            nrows -= 1
            data = pd.read_csv(path, encoding='iso-8859-1', sep="\t", skiprows=nrows)
        data['ptid'] = path.split("/")[-1].split("-")[-1].split("_")[0]
        # print(f"path: {path.split('/')[-1]}, ptid: {data.ptid.iloc[0]}")
        data.ptid = data.ptid.astype(int)
        data['filter_pass'] = path.split("/")[-1].split(".")[0].split("_")[-1].strip()
        data['filename'] = path.split("/")[-1]
        usecols = ['ptid', 'Name', 'Block', 'filter_pass', 'F532 Mean', 'B532 Mean', 'filename']
        data = data.loc[~data.Name.isin(["empty","Empty","DyeSpot"]),usecols]
        data = data.rename(columns={'Name':'spot_name',
                                    'F532 Mean':'f532_mean',
                                    'B532 Mean':'b532_mean',
                                    'Block':'block'})
        return data

    data = pd.DataFrame()
    for file in files:
        data = pd.concat([data, get_data(datadir + file)])

    # from the lab: we only want the "Hi" data
    data = data.loc[data.filter_pass=="Hi"]

    metadata_path = "/trials/vaccine/p302/s001/qdata/LabData/AE_assays_pass-through/Glycan_array/IgE/lab_processed_data/HVTN302-PatientSample-GlycanAnalysis.xlsx"
    metadata = pd.read_excel(metadata_path)

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

    # drop rows for which no guspec -- these wouldve been blanks on slides
    # software reads it anyway
    data = data.loc[data.guspec.notna()]

    ## read glycan info ------------------------------------------------------##
    glycan_info_path = "/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN302/assays/AE_assays/glycan_microarray/misc_files/Scheif-GlycanArrayList-Final.xlsx"
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
    data['f532_mean_minus_b532_mean'] = data.f532_mean - data.b532_mean

    ## read in ldms ----------------------------------------------------------##
    ldms = pd.read_csv(constants.LDMS_PATH_HVTN, usecols=constants.STANDARD_COLS, dtype=constants.LDMS_DTYPE_MAP)
    ldms = ldms.loc[ldms.lstudy==302.]
    ldms = ldms.loc[ldms.guspec.isin(data.guspec)]

    ## hand-entered data -----------------------------------------------------##
    metadata_dict = {
        "network": "HVTN",
        "upload_lab_id": "P4",
        "assay_lab_name": "Paulson (Scripps)",
        "assay_type": "Glycan Microarray",
        "assay_precision": "Semi-Quantitative",
        "assay_details":"Custom glycan layout",
        "instrument": "Innoscan 1100AL",
        "lab_software_version": " Mapix",
        "specrole": "Sample",
        "glycan_mapping_file": "Scheif-GlycanArrayList-Final.xlsx",
        "sample_mapping_file": "HVTN302-PatientSample-GlycanAnalysis.xlsx"
    }

    ## standard processing ---------------------------------------------------##
    outputs = sdmc.standard_processing(
        input_data=data.drop(columns='ptid'),
        input_data_path=datadir + files[0],
        guspec_col="guspec",
        network="HVTN",
        metadata_dict=metadata_dict,
        ldms=ldms
    )

    # standard_processing not handling variable input file name
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
        'f532_mean',
        'b532_mean',
        'f532_mean_minus_b532_mean',
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

    ## save results ----------------------------------------------------------##
    today = datetime.date.today().isoformat()
    fname = f"DRAFT_HVTN302_IgE_glycan_data_processed_{today}.txt"
    savedir = "/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN302/assays/AE_assays/glycan_microarray/misc_files/data_processing/IgE"
    outputs.to_csv(savedir + fname, sep="\t", index=False)

    # manifest comparison ----------------------------------------------------##
    manifest_path = "/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN302/assays/AE_assays/glycan_microarray/misc_files/hvtn302_glycan_shipping_manifest_from_nick.txt"
    manifest = pd.read_csv(manifest_path, sep="\t")
    mismatch = set(outputs.guspec).symmetric_difference(manifest.GLOBAL_ID)
    if len(mismatch) > 0:
        raise Exception("guspecs don't match manifest")

def pivot():
    # pivot summaries --------------------------------------------------------##
    guspec_glycan = pd.pivot_table(outputs, index='guspec', columns='glycan_m_number', aggfunc='count')['f532_mean_minus_b532_mean']
    ptid_visitno_glycan = pd.pivot_table(outputs, index=['ptid','visitno'], columns='glycan_m_number', aggfunc='count')['f532_mean_minus_b532_mean']
    ptid_visitno = pd.pivot_table(outputs, index='ptid', columns='visitno', fill_value=0, aggfunc='count')['f532_mean_minus_b532_mean']

    guspec_glycan.to_excel(savedir + "HVTN302_IgE_glycan_microarray_guspec_glycan_summary.xlsx")
    ptid_visitno_glycan.to_excel(savedir + "HVTN302_IgE_glycan_microarray_ptid_visitno_glycan_summary.xlsx")
    ptid_visitno.to_excel(savedir + "HVTN302_IgE_glycan_microarray_ptid_visitno_summary.xlsx")

if __name__=="__main__":
    main()
    # pivot()
