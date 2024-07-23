## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 07/15/2024
# Purpose: For CAVD VISC: Concatenate / merge Glycan microarray data
## ---------------------------------------------------------------------------##
import pandas as pd
import numpy as np
import os
import datetime
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
    input_data = pd.DataFrame()
    for site in os.listdir(parentdir):
        datadir = parentdir + "/"  + site + "/Processed Data/"
        files = os.listdir(datadir)
        files = [f for f in files if f[-3:]=="txt"]
        sub = pd.DataFrame()
        for file in files:
            sub = pd.concat([sub, get_data(datadir + file)])
        sub["site"] = site
        input_data = pd.concat([input_data, sub])

    # verify that there are 6 rows per unique combo
    t = input_data.groupby(['ptid','site','spot_name','filter_pass','block']).count().reset_index()
    count_per = t.f488_mean.unique().tolist()
    if len(count_per) > 1 or count_per[0]!=6:
        print(f"t.f488_mean.unique(): {t.f488_mean.unique()}")
        raise Exception("There should be 6 rows per unique spit/site/spot name/filter/block!")

    # subset to 'hi' pass only, per the lab
    input_data = input_data.loc[input_data.filter_pass=="Hi"]

    # add primary result column
    data = input_data.copy()
    data['f488_mean_minus_b488_mean'] = np.round(input_data.f488_mean - input_data.b488_mean, 2)

    # Merge on timepoint column ----------------------------------------------##
    sample_id_path = '/networks/vtn/lab/SDMC_labscience/operations/documents/templates/assay/template_testing/IAVI-G002_ArrayDataKey.xlsx'
    sample_id_df = pd.read_excel(sample_id_path)

    # make sure ptid/visit uniquely identifies
    if len(sample_id_df[['PTID','Visit']].drop_duplicates())!=len(sample_id_df):
        raise Exception("ptid/visit arent a unique identifier in lab sample mapping!")

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

    # Merge on CAVD metadata (incl sample identifier) ----------------------------------------------##
    cavd_sample_list_path = "/networks/cavd/VDCs/Schief/Schief_856-G002/SkinReactions/samples/Skin AE sample testing.xlsx"
    cavd_sample_list = pd.read_excel(cavd_sample_list_path)
    cavd_sample_list = cavd_sample_list.loc[cavd_sample_list.Destination=="GCF"]

    data.Visit = data.Visit.astype(int)

    cavd_sample_list = cavd_sample_list.rename(columns = {"Subject ID": "ptid", "Study Day":"study_week"})
    cavd_sample_list['Date Drawn'] = cavd_sample_list['Date Drawn'].dt.date

    data = data.merge(cavd_sample_list, on=["Study ID", "ptid", "Visit", "study_week"], how="left")

    # verify that there are still 6 rows per unique combo
    t = data.groupby(['Current Label', 'spot_name', 'block']).count().reset_index()
    count_per = t.f488_mean.unique().tolist()
    if len(count_per) > 1 or count_per[0]!=6:
        print(f"t.f488_mean.unique(): {t.f488_mean.unique()}")
        raise Exception("There should be 6 rows per unique sample/spot_name/block!")

    compare = cavd_sample_list[['ptid','study_week']].merge(
        sample_id_df,
        left_on=['ptid', 'study_week'],
        right_on=['PTID',"Study Week"], how = 'outer'
    )

    if len(compare)!=len(cavd_sample_list):
        raise Exception("lab sample list and cavd sample list don't match!")

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
        "lab_sample_mapping_file": ["IAVI-G002_ArrayDataKey.xlsx"],
        "cavd_sample_mapping_file": ["Skin AE sample testing.xlsx"],
        "fred_hutch_processing_timestamp": [fred_hutch_processing_timestamp],
        "network": "CAVD",
    }

    data = data.merge(pd.DataFrame(metadata_dict), how='cross')

    ## standard formatting ---------------------------------------------------##
    data.columns = [i.lower().replace(" ","_") for i in data.columns]

    data = data.rename(columns = {
        'filename':'input_file_name',
        'blocks': 'block_group',
        'study_id': 'lab_submitted_study_id',
        'current_label': 'sample_id',
        'volume': 'sample_volume_total',
        'volume_unit': 'sample_volume_unit'
    })

    reorder = [
        'network',
        'lab_submitted_study_id',
        'sample_id',
        'ptid',
        'visit',
        'study_week',
        'date_drawn',
        'block',
        'block_group',
        'material_type',
        'material_modifiers',
        'specrole',
        'upload_lab_id',
        'assay_lab_name',
        'assay_type',
        'instrument',
        'lab_software_version',
        'assay_precision',
        'assay_details',
        'site',
        'group',
        'destination',
        'requisition_id',
        'spot_name',
        'glycan_m_number',
        'glycan_structure',
        'f488_mean',
        'b488_mean',
        'f488_mean_minus_b488_mean',
        'filter_pass',
        'sample_volume_total',
        'sample_volume_unit',
        'input_file_name',
        'glycan_mapping_file',
        'lab_sample_mapping_file',
        'cavd_sample_mapping_file',
        'input_data_timestamp',
        'fred_hutch_processing_timestamp'
    ]

    mismatch = set(reorder).symmetric_difference(data.columns)
    if len(mismatch) > 0:
        raise Exception(f"Trying to drop or add columns: {mismatch}")

    outputs = data[reorder]

    outputs.sample_id = outputs.sample_id.astype(str)

    ## save results ----------------------------------------------------------##
    savedir = "/networks/cavd/VDCs/Schief/Schief_856-G002/SkinReactions/data/Glycan_array_Scripps/processed_data/"
    today = datetime.date.today().isoformat()
    outputs.to_csv(savedir + f"DRAFT_CAVD_G002_Glycan_Microarray_data_processed_{today}.txt", sep="\t", index=False)

if __name__=="__main__":
    main()
