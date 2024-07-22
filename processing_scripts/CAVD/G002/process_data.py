## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 07/15/2024
# Purpose: For CAVD VISC: Concatenate / merge Glycan microarray data
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
    def get_data(path):
        nrows = 35
        data = pd.read_csv(path, encoding='iso-8859-1', sep="\t", skiprows=nrows)
        while 'Block' not in data.columns and nrows > 0:
            nrows -= 1
            data = pd.read_csv(path, encoding='iso-8859-1', sep="\t", skiprows=nrows)
        data['ptid'] = path.split("/")[-1].split("_")[0]
        # data.ptid = data.ptid.astype(int)
        data['filter_pass'] = path.split("/")[-1].split(".")[0].split("_")[-1]
        usecols = ['ptid', 'Name', 'Block', 'filter_pass', 'F488 Mean', 'B488 Mean']
        data = data.loc[~data.Name.isin(["empty","Empty","DyeSpot"]),usecols]
        data = data.rename(columns={'Name':'spot_name',
                                    'F488 Mean':'f488_mean',
                                    'B488 Mean':'b488_mean',
                                    'Block':'block'})
        input_data_timestamp = datetime.datetime.fromtimestamp(os.path.getmtime(path)).replace(microsecond=0).isoformat()
        data['input_data_timestamp'] = input_data_timestamp
        return data

    parentdir = '/networks/cavd/VDCs/Schief/Schief_856-G002/SkinReactions/data/Glycan_array_Scripps/20240327-01/OneDrive_1_3-27-2024/'
    data = pd.DataFrame()
    for site in os.listdir(parentdir):
        datadir = parentdir + "/"  + site + "/Processed Data/"
        files = os.listdir(datadir)
        files = [f for f in files if f[-3:]=="txt"]
        sub = pd.DataFrame()
        for file in files:
            sub = pd.concat([sub, get_data(datadir + file)])
        sub["site"] = site
        data = pd.concat([data, sub])

    # verify that there are 6 rows per unique combo
    t = data.groupby(['ptid','site','spot_name','filter_pass','block']).count().reset_index()
    if not t.f488_mean.unique().tolist==[6]:
        raise Exception("There should be 6 rows per unique sample/spot name/filter/block!")

    # add primary result column
    data['f488_mean_minus_b488_mean'] = data.f488_mean - data.b488_mean

    # subset to 'hi' pass only, per the lab
    data = data.loc[data.filter_pass=="Hi"]

    # Merge on timepoint column ----------------------------------------------##
    sample_id_path = '/networks/vtn/lab/SDMC_labscience/operations/documents/templates/assay/template_testing/IAVI-G002_ArrayDataKey.xlsx'
    sample_id_df = pd.read_excel(sample_id_path)

    # verify that these are an exact match
    if len(set(sample_id_df.PTID).symmetric_difference(data.ptid)) > 0:
        raise Exception("Sample list PTIDS don't match ptids in data")

    # make sure none of the sample mapping rows are NaN
    if len(sample_id_df.loc[sample_id_df.PTID.isna() | sample_id_df['Study Week'].isna() | sample_id_df.Filename.isna()]) > 0:
        raise Exception("NaNs in sample mapping file")

    # map sample info onto data using ptid/timepoint
    block_to_timepoint = {**{i:"Wk 0" for i in range(13)}, **{i:"Wk 8" for i in range(13,25)}, **{i:"Wk 10" for i in range(25,37)}}
    data['study_week'] = data.block.map(block_to_timepoint)
    data = data.merge(sample_id_df, left_on=['ptid','study_week'], right_on=['PTID','Study Week'], how='left')

    # note the ptids for which we don't have three timepoints
    missings = data.loc[data.Filename.isna(),['ptid','study_week']].drop_duplicates()
    if len(missings) > 0:
        print("The following ptid/study week combos are missing from the sample list from the lab:")
        print(missings)

    # dropping all rows that don't correspond to something in the sample manifest
    data = data.loc[~data.Filename.isna()]

    # make sure 'filename' row is consistent
    if len(data.loc[data.ptid!=data.Filename.str[:-9]]) > 0:
        raise Exception("Ptid mismatch identified")

    # site and Origin should be a match:
    print("site and Origin should be a match:")
    print(data[['site','Origin']].drop_duplicates())

    # drop duplicate columns
    data = data.drop(columns = ['Origin', 'PTID', 'Study Week'])

    ## merge on glycan info --------------------------------------------------##
    glycan_info_path = "/networks/cavd/VDCs/Schief/Schief_856-G002/SkinReactions/data/Glycan_array_Scripps/Scheif-GlycanArrayList-Final.xlsx"
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

    ## merge on metadata -----------------------------------------------------##
    fred_hutch_processing_timestamp = datetime.datetime.now().replace(microsecond=0).isoformat()

    metadata_dict = {
        "network": ["HVTN"],
        "upload_lab_id": ["P4"],
        "assay_lab_name": ["Paulson (Scripps)"],
        "assay_type": ["Glycan Microarray"],
        "assay_precision": ["Semi-Quantitative"],
        "assay_details":["Custom glycan layout"],
        "instrument": ["Innoscan 1100AL"],
        "lab_software_version": ["Mapix"],
        "specrole": ["Sample"],
        "glycan_mapping_file": ["Scheif-GlycanArrayList-Final.xlsx"],
        "sample_mapping_file": ["IAVI-G002_ArrayDataKey.xlsx"],
        "fred_hutch_processing_timestamp": [fred_hutch_processing_timestamp],
        "network": "CAVD",
        "sample_id": "NA",
    }

    data = data.merge(pd.DataFrame(metadata_dict), how='cross')

    ## standard formatting ---------------------------------------------------##
    data.columns = [i.lower().replace(" ","_") for i in data.columns]

    data = data.rename(columns = {
        'filename':'input_file_name',
        'blocks': 'block_group',
        'study_id': 'lab_submitted_study_id'
    })

    reorder = [
        'network',
        'lab_submitted_study_id',
        'sample_id',
        'ptid',
        'visit',
        'study_week',
        'block_group',
        'block',
        'specrole',
        'upload_lab_id',
        'assay_lab_name',
        'assay_type',
        'instrument',
        'lab_software_version',
        'assay_precision',
        'assay_details',
        'site',
        'spot_name',
        'glycan_m_number',
        'glycan_structure',
        'f488_mean',
        'b488_mean',
        'f488_mean_minus_b488_mean',
        'filter_pass',
        'input_file_name',
        'glycan_mapping_file',
        'sample_mapping_file',
        'input_data_timestamp',
        'fred_hutch_processing_timestamp'
    ]

    mismatch = set(reorder).symmetric_difference(data.columns)
    if len(mismatch) > 0:
        raise Exception(f"Trying to drop or add columns: {mismatch}")

    outputs = data[reorder]

    ## save results ----------------------------------------------------------##
    savedir = "/networks/vtn/lab/SDMC_labscience/studies/VISC/G002/data_processing/"
    timestamp = datetime.date.today().isoformat()
    outputs.to_csv(savedir + "EXAMPLE_CAVD_G002_Glycan_Microarray_data_processed_{today}.txt", sep="\t", index=False)

if __name__=="__main__":
    main()
