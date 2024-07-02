## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 06/18/2024
# Purpose:  - Concatenate / merge Glycan microarray data for stats
## ---------------------------------------------------------------------------##
import pandas as pd
import numpy as np
import os
import datetime
import yaml
import sdmc_tools.process as sdmc
import sdmc_tools.constants as constants
## ---------------------------------------------------------------------------##

def main():
    ## read in raw data ------------------------------------------------------##
    datadir = '/trials/vaccine/p302/s001/qdata/LabData/AE_assays_pass-through/Glycan_array/lab_processed_data/'
    files = os.listdir(datadir)
    files = [f for f in files if f[-3:]=="txt"]

    # read in data
    def get_data(path):
        data = pd.read_csv(path, encoding='iso-8859-1', sep="\t", skiprows=35)
        data['ptid'] = path.split("/")[-1].split("_")[0]
        data.ptid = data.ptid.astype(int)
        data['filter_pass'] = path.split("/")[-1].split(".")[0].split("_")[-1]
        usecols = ['ptid', 'Name', 'Block', 'filter_pass', 'F488 Mean', 'B488 Mean']
        data = data.loc[~data.Name.isin(["empty","Empty","DyeSpot"]),usecols]
        return data

    data = pd.DataFrame()
    for file in files:
        data = pd.concat([data, get_data(datadir + file)])

    # from the lab: we only want the "Hi" data
    data = data.loc[data.filter_pass=="Hi"]

    ## read in metadata ------------------------------------------------------##
    metadata_path = "/trials/vaccine/p302/s001/qdata/LabData/AE_assays_pass-through/Glycan_array/lab_processed_data/metadata/HVTN302-PatientSample-GlycanAnalysis.xlsx"
    metadata = pd.read_excel(metadata_path)

    metadata.columns = [i.lower().replace(" ","_") for i in metadata.columns]
    metadata = metadata.rename(columns = {'original_id':'guspec', 'subject_id':'ptid', 'box_':'box'})
    metadata['visitno'] = metadata['patient-visit'].str.split("-V", expand=True)[1].astype(int)

    ## add 'timepoint' to both data and metadata, then merge -----------------##
    def label_timepoint(ptid, visitno):
        visits = metadata.loc[(metadata.ptid==ptid)].visitno.tolist()
        return visits.index(visitno)

    metadata['timepoint'] = metadata.apply(lambda x: label_timepoint(x['ptid'], x['visitno']), axis=1)

    block_to_timepoint = {**{i:0 for i in range(13)}, **{i:1 for i in range(13,25)}, **{i:2 for i in range(25,37)}}
    data["timepoint"] = data.Block.map(block_to_timepoint)
    data = data.merge(metadata[['ptid','timepoint','guspec']], on=['ptid','timepoint'], how='left')

    # dropping merge column
    data = data.drop(columns='timepoint')

    # dropping 8 rows for which don't have data. it looks like we never sent those samples
    data = data.loc[data.guspec.notna()]

    ## read glycan info ------------------------------------------------------##
    glycan_info_path = "/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN302/assays/AE_assays/glycan_microarray/misc_files/Scheif-GlycanArrayList-Final.xlsx"
    glyc_info = pd.read_excel(glycan_info_path)
    glyc_info = glyc_info.rename(columns={'Sample':'s', 'M-Number':'ref_glycan'})

    # merge on M-Number from glyc_info
    tmp = data[['Name']].drop_duplicates()
    tmp['s'] = tmp.Name.str.partition("-", expand=True)[0]
    tmp.s = tmp.s.astype(int)
    tmp['data_glycan'] = tmp.Name.str.partition("-", expand=True)[2]

    tmp = glyc_info.merge(tmp, on='s', how='outer')
    Name_to_glycan_m_number_map = {data:ref for (data,ref) in zip(tmp.Name, tmp.ref_glycan)}
    data['m_number'] = data.Name.map(Name_to_glycan_m_number_map)
    data.loc[data.Name=="038-418Sp19", "Name"] = "038-419Sp19"

    ## final pre-processing --------------------------------------------------##
    # standardize columns
    data.columns = [i.lower().replace(" ","_") for i in data.columns]

    # add primary result column
    data['f488_mean_minus_b488_mean'] = data.f488_mean - data.b488_mean

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
        "m_number_source": "Scheif-GlycanArrayList-Final.xlsx",
        "guspec_source": "HVTN302-PatientSample-GlycanAnalysis.xlsx"
    }

    ## standard processing ---------------------------------------------------##
    outputs = sdmc.standard_processing(
        input_data=data.drop(columns='ptid'),
        input_data_path=datadir + "604369981_G_Hi.txt",
        guspec_col="guspec",
        network="HVTN",
        metadata_dict=metadata_dict,
        ldms=ldms
    )

    # input file is specific to each row; need to fix
    outputs['input_file_name'] = outputs.ptid + "_G_Hi.txt"

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
        'name',
        'm_number',
        'f488_mean',
        'b488_mean',
        'f488_mean_minus_b488_mean',
        'filter_pass',
        'block',
        'guspec_source',
        'm_number_source',
        'input_file_name',
        'sdmc_processing_datetime',
        'sdmc_data_receipt_datetime',
    ]

    outputs = outputs[reorder]

    ## save results ----------------------------------------------------------##
    today = datetime.date.today().isoformat()
    fname = f"DRAFT_HVTN302_glycan_data_processed_{today}.txt"
    savedir = "/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN302/assays/AE_assays/glycan_microarray/misc_files/data_processing/"

    outputs.to_csv(savedir + fname, sep="\t", index=False)

def pivot():
    savedir = "/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN302/assays/AE_assays/glycan_microarray/misc_files/data_processing/"
    outputs = pd.read_csv(savedir + "HVTN302_glycan_data_processed_2024-06-18.txt", sep="\t")

    guspec_by_glycan = pd.pivot_table(
                            outputs,
                            index='guspec',
                            columns='glycan_m_number',
                            values='f488_mean_minus_b488_mean',
                            aggfunc='count'
                            )
    guspec_by_glycan.to_excel(savedir + "HVTN302_Glycan_Microarray_guspec_x_glycan_summary_count.xlsx")

    ptid_visit_glycan = pd.pivot_table(
                            outputs,
                            index=['ptid','visitno'],
                            columns='glycan_m_number',
                            aggfunc='count'
                            )
    ptid_visit_glycan.to_excel(savedir + "HVTN302_Glycan_Microarray_ptid_visitno_x_glycan_summary_count.xlsx")

if __name__ == '__main__':
    main()
    pivot()
