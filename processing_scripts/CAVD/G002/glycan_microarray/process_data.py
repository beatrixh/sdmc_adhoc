## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 10/16/2024
# Revision: 01/21/2025
# Purpose: For CAVD VISC: Concatenate / merge Glycan microarray data
## ---------------------------------------------------------------------------##
import pandas as pd
import numpy as np
import os
import datetime

def main():
    ## read in and concat results data ---------------------------------------##
    def get_data(path, channel):
        nrows = 35
        data = pd.read_csv(path, encoding='iso-8859-1', sep="\t", skiprows=nrows)
        while 'Block' not in data.columns and nrows > 0:
            nrows -= 1
            data = pd.read_csv(path, encoding='iso-8859-1', sep="\t", skiprows=nrows)
        data['ptid'] = path.split("/")[-1].split("_")[0]
        data['filename'] = path.split("/")[-1]
        data['filter_pass'] = path.split("/")[-1].split(".")[0].split("_")[-1]
        data['signal_channel'] = channel
        input_data_timestamp = datetime.datetime.fromtimestamp(os.path.getmtime(path)).replace(microsecond=0).isoformat()
        data['input_data_timestamp'] = input_data_timestamp

        data = data.rename(columns={'Name':'spot_name',
                                    f'F{channel} Mean': f'mean_signal',
                                    f'B{channel} Mean': f'mean_background_signal',
                                    'Block':'block'})
        usecols = ['ptid', 'spot_name', 'block', 'filter_pass', 'signal_channel', 'mean_signal', 'mean_background_signal', 'filename', 'input_data_timestamp']

        data = data.loc[~data.spot_name.isin(["empty","Empty","DyeSpot"]),usecols]
        return data

    ## IgM ##
    parentdir = "/networks/cavd/VDCs/Schief/Schief_856-G002/SkinReactions/data/Glycan_array_Scripps/20241011-01/"
    igm = pd.DataFrame()
    for site in os.listdir(parentdir):
        if os.path.isdir(parentdir + site):
            datadir = parentdir + site + "/ProcessedData/"
            files = os.listdir(datadir)
            files = [f for f in files if f[-3:]=="txt"]
            sub = pd.DataFrame()
            for file in files:
                sub = pd.concat([sub, get_data(datadir + file, 532)])
            sub["site"] = site
            igm = pd.concat([igm, sub])
    igm = igm.loc[igm.filter_pass=="Lo"]
    igm['isotype'] = 'IgM'
    igm['background_subtraced_mean_signal'] = igm['mean_signal'] - igm['mean_background_signal']
    # floating point error
    igm['background_subtraced_mean_signal'] = igm['background_subtraced_mean_signal'].round(2)
    igm.site = igm.site.map({'Emory': 'Emory', 'FH': 'Fred Hutch', 'GW': 'George Washington', 'UT': 'University of Texas'})

    ## IgE ##
    parentdir = "/networks/cavd/VDCs/Schief/Schief_856-G002/SkinReactions/data/Glycan_array_Scripps/20240903-01/"
    ige = pd.DataFrame()
    for site in os.listdir(parentdir):
        if os.path.isdir(parentdir + site):
            datadir = parentdir + site + "/ProcessedData/"
            files = os.listdir(datadir)
            files = [f for f in files if f[-3:]=="txt"]
            sub = pd.DataFrame()
            for file in files:
                sub = pd.concat([sub, get_data(datadir + file, 532)])
            sub["site"] = site
            ige = pd.concat([ige, sub])

    ige = ige.loc[ige.filter_pass=="Hi"]
    ige['isotype'] = 'IgE'
    ige['background_subtraced_mean_signal'] = ige['mean_signal'] - ige['mean_background_signal']
    # floating point error
    ige['background_subtraced_mean_signal'] = ige['background_subtraced_mean_signal'].round(2)
    ige.site = ige.site.map({'Emory': 'Emory', 'FH': 'Fred Hutch', 'GW': 'George Washington', 'UT': 'University of Texas'})

    ## IgG ##
    parentdir = '/networks/cavd/VDCs/Schief/Schief_856-G002/SkinReactions/data/Glycan_array_Scripps/20240327-01/OneDrive_1_3-27-2024/'
    igg = pd.DataFrame()
    for site in os.listdir(parentdir):
        if os.path.isdir(parentdir + site):
            datadir = parentdir  + site + "/Processed Data/"
            files = os.listdir(datadir)
            files = [f for f in files if f[-3:]=="txt"]
            sub = pd.DataFrame()
            for file in files:
                sub = pd.concat([sub, get_data(datadir + file, 488)])
            sub["site"] = site
            igg = pd.concat([igg, sub])

    igg = igg.loc[igg.filter_pass=="Hi"]
    igg['isotype'] = "IgG"
    igg['background_subtraced_mean_signal'] = igg['mean_signal'] - igg['mean_background_signal']
    # floating point error
    igg['background_subtraced_mean_signal'] = igg['background_subtraced_mean_signal'].round(2)

    ## Merge timepoint data (lab sample list) --------------------------------##
    def merge_on_visitno(data, visitno_df):
        # make sure ptid/visit uniquely identifies
        if len(visitno_df[['PTID','Visit']].drop_duplicates())!=len(visitno_df):
            raise Exception("ptid/visit arent a unique identifier in lab sample mapping!")

        # verify that these are an exact match
        if len(set(visitno_df.PTID).symmetric_difference(data.ptid)) > 0:
            raise Exception("Sample list PTIDS don't match ptids in data")

        # make sure none of the sample mapping rows are NaN
        if len(visitno_df.loc[visitno_df.PTID.isna() | visitno_df['Study Week'].isna() | visitno_df.Filename.isna()]) > 0:
            raise Exception("NaNs in sample mapping file")

        # make sure 'filename' row is consistent
        if len(visitno_df.loc[visitno_df.PTID!=visitno_df.Filename.str[:-9]]) > 0:
            raise Exception("Ptid mismatch identified")
        visitno_df = visitno_df.drop(columns='Filename')

        # map sample info onto data using ptid/timepoint
        block_to_timepoint = {**{i:"Wk 0" for i in range(13)}, **{i:"Wk 8" for i in range(13,25)}, **{i:"Wk 10" for i in range(25,37)}}
        data['study_week'] = data.block.map(block_to_timepoint)
        data = data.merge(visitno_df, left_on=['ptid','study_week'], right_on=['PTID','Study Week'], how='left')

        # note the ptids for which we don't have three timepoints
        missings = data.loc[data.Visit.isna(),['ptid','study_week']].drop_duplicates()
        if len(missings) > 0:
            print("The following ptid/study week combos are missing from the sample list from the lab:")
            print(missings)

        # dropping all rows that don't correspond to something in the sample manifest
        data = data.loc[~data.Visit.isna()]

        # site and Origin should be a match:
        print("site and Origin should be a match:")
        print(data[['site','Origin']].drop_duplicates())

        return data

    ## IgM ##
    igm_visitno_metadata_path = '/networks/cavd/VDCs/Schief/Schief_856-G002/SkinReactions/data/Glycan_array_Scripps/20241011-01/IAVI-G002_ArrayDataKey.xlsx'
    igm_visitno_df = pd.read_excel(igm_visitno_metadata_path)
    igm_visitno_df = igm_visitno_df.drop(columns=[
        'Filename (IgG)',
        'Blocks',
        'Filename (IgE)',
        'Blocks.1',
        'Unnamed: 14'
    ])
    igm_visitno_df = igm_visitno_df.dropna(axis=1)
    igm_visitno_df = igm_visitno_df.rename(columns={
        'Filename (IgM)': 'Filename',
        'Blocks.2': 'Blocks'
    })

    igm2 = merge_on_visitno(igm, igm_visitno_df)
    igm2 = igm2.drop(columns = ['Origin', 'PTID', 'Study Week'])

    ## IgG ##
    igg_visitno_metadata_path = '/networks/cavd/VDCs/Schief/Schief_856-G002/SkinReactions/data/Glycan_array_Scripps/IAVI-G002_ArrayDataKey_19Jul2024.xlsx'
    igg_visitno_df = pd.read_excel(igg_visitno_metadata_path)
    igg2 = merge_on_visitno(igg, igg_visitno_df)
    igg2 = igg2.drop(columns = ['Origin', 'PTID', 'Study Week'])

    ## IgE ##
    ige_visitno_metadata_path = '/networks/cavd/VDCs/Schief/Schief_856-G002/SkinReactions/data/Glycan_array_Scripps/IAVI-G002_ArrayDataKey_10Sep2024.xlsx'
    ige_visitno_df = pd.read_excel(ige_visitno_metadata_path)
    ige_visitno_df = ige_visitno_df.dropna(axis=1)
    ige_visitno_df = ige_visitno_df.drop(columns=['Filename (IgG)', 'Blocks'])
    ige_visitno_df = ige_visitno_df.rename(columns={
        'Filename (IgE)': 'Filename',
        'Blocks.1': 'Blocks'
    })
    ige2 = merge_on_visitno(ige, ige_visitno_df)
    ige2 = ige2.drop(columns = ['Origin', 'PTID', 'Study Week'])

    ## Merge sammple id (cavd sample list) -----------------------------------##
    cavd_sample_list_path = "/networks/cavd/VDCs/Schief/Schief_856-G002/SkinReactions/samples/Skin AE sample testing.xlsx"
    def merge_on_sample_identifier(data, cavd_sample_list_path):
        cavd_sample_list = pd.read_excel(cavd_sample_list_path)
        cavd_sample_list = cavd_sample_list.loc[cavd_sample_list.Destination=="GCF"]

        data.Visit = data.Visit.astype(int)

        cavd_sample_list = cavd_sample_list.drop(columns=["Volume", "Volume Unit"])
        cavd_sample_list = cavd_sample_list.rename(columns = {"Subject ID": "ptid", "Study Day":"study_week"})
        cavd_sample_list['Date Drawn'] = cavd_sample_list['Date Drawn'].dt.date

        data = data.merge(cavd_sample_list, on=["Study ID", "ptid", "Visit", "study_week"], how="left")

        # verify that there are still 6 rows per unique combo
        t = data.groupby(['Current Label', 'spot_name', 'block']).count().reset_index()
        count_per = t.mean_signal.unique().tolist()
        if len(count_per) > 1 or count_per[0]!=6:
            print(f"t.mean_signal.unique(): {t.mean_signal.unique()}")
            raise Exception("There should be 6 rows per unique sample/spot_name/block!")
        return data

    ## IgM ##
    igm3 = merge_on_sample_identifier(igm2, cavd_sample_list_path)

    ## IgG ##
    igg3 = merge_on_sample_identifier(igg2, cavd_sample_list_path)

    ## IgE ##
    ige3 = merge_on_sample_identifier(ige2, cavd_sample_list_path)

    ## make sure metadata matches one another --------------------------------##
    def compare_metadata(cavd_sample_list_path, visitno_metadata_path):
        cavd_sample_list = pd.read_excel(cavd_sample_list_path)
        visitno_df = pd.read_excel(visitno_metadata_path)

        cavd_sample_list['Subject ID'] = cavd_sample_list['Subject ID'].astype(str)
        visitno_df.PTID = visitno_df.PTID.astype(str)

        compare = cavd_sample_list[['Subject ID','Study Day']].merge(
            visitno_df,
            left_on=['Subject ID','Study Day'],
            right_on=['PTID', 'Study Week'], how = 'outer'
        )

        if len(compare)!=len(cavd_sample_list):
            raise Exception("lab sample list and cavd sample list don't match!")

    compare_metadata(cavd_sample_list_path, igg_visitno_metadata_path)
    compare_metadata(cavd_sample_list_path, ige_visitno_metadata_path)
    compare_metadata(cavd_sample_list_path, igm_visitno_metadata_path)

    ## merge on glycan metadata ----------------------------------------------##
    glycan_info_path = "/networks/cavd/VDCs/Schief/Schief_856-G002/SkinReactions/data/Glycan_array_Scripps/Scheif-GlycanArrayList-Final.xlsx"

    def merge_on_glycan_info(data, glycan_info_path):
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
        data.columns = [i.lower().replace(" ","_") for i in data.columns]
        return data

    ## IgM ##
    igm4 = merge_on_glycan_info(igm3, glycan_info_path)

    ## IgG ##
    igg4 = merge_on_glycan_info(igg3, glycan_info_path)

    ## IgE ##
    ige4 = merge_on_glycan_info(ige3, glycan_info_path)

    ## standard processing ---------------------------------------------------##
    ## merge on metadata ##
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
        "lab_sample_mapping_file": [igm_visitno_metadata_path.split("/")[-1]],
        "cavd_sample_mapping_file": ["Skin AE sample testing.xlsx"],
        "fred_hutch_processing_timestamp": [fred_hutch_processing_timestamp],
        "network": "CAVD",
    }

    igm5 = igm4.merge(pd.DataFrame(metadata_dict), how='cross')

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
        "lab_sample_mapping_file": [igg_visitno_metadata_path.split("/")[-1]],
        "cavd_sample_mapping_file": ["Skin AE sample testing.xlsx"],
        "fred_hutch_processing_timestamp": [fred_hutch_processing_timestamp],
        "network": "CAVD",
    }

    igg5 = igg4.merge(pd.DataFrame(metadata_dict), how='cross')

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
        "lab_sample_mapping_file": [ige_visitno_metadata_path.split("/")[-1]],
        "cavd_sample_mapping_file": ["Skin AE sample testing.xlsx"],
        "fred_hutch_processing_timestamp": [fred_hutch_processing_timestamp],
        "network": "CAVD",
    }
    ige5 = ige4.merge(pd.DataFrame(metadata_dict), how='cross')

    ## standard formatting ##
    data = pd.concat([igg5, ige5, igm5])

    data.columns = [i.lower().replace(" ","_") for i in data.columns]

    data = data.rename(columns = {
        'filename':'input_file_name',
        'blocks': 'block_group',
        'study_id': 'lab_submitted_study_id',
        'current_label': 'sample_id',
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
        'isotype',
        'spot_name',
        'glycan_m_number',
        'glycan_structure',
        'signal_channel',
        'mean_signal',
        'mean_background_signal',
        'background_subtraced_mean_signal', ## REACH OUT TO STATS TO UPDATE
        'filter_pass',
        'input_file_name',
        'glycan_mapping_file',
        'lab_sample_mapping_file',
        'cavd_sample_mapping_file',
        'input_data_timestamp',
        'fred_hutch_processing_timestamp',
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

    pivot_summary = pd.pivot_table(outputs, index='ptid',columns=['isotype','visit'], aggfunc='count', fill_value=0)['mean_signal']
    pivot_summary.to_excel(savedir + f"CAVD_G002_glycan_ptid_visit_isotype_summary_{today}.xlsx")

    ## check against 10/2024 version -----------------------------------------##

    # IgM signal_channel, laser power ('filter_pass') and results have changed
    # processing datetime has changed
    # everything else should be the same

    old = pd.read_csv(savedir + 'DRAFT_CAVD_G002_Glycan_Microarray_data_processed_2024-10-16.txt', sep="\t")

    # correct floating point error
    old.background_subtraced_mean_signal = old.background_subtraced_mean_signal.round(2)

    # pandas types reading in differently; correct
    outputs.sample_id = outputs.sample_id.astype(str)
    old.sample_id = old.sample_id.astype(str)

    outputs.date_drawn = outputs.date_drawn.astype(str)
    old.date_drawn = old.date_drawn.astype(str)

    outputs.signal_channel = outputs.signal_channel.astype(float)
    old.signal_channel = old.signal_channel.astype(float)

    # align to comparable format
    old = old[list(outputs.columns)]

    outputs = outputs.sort_values(by=list(outputs.columns)).reset_index(drop=True)
    old = old.sort_values(by=list(outputs.columns)).reset_index(drop=True)

    ## compare:
    # after correcting for floating point math error
    # IgG and IgE are identical!
    outputs.loc[outputs.isotype!="IgM"].reset_index(drop=True).compare(old.loc[old.isotype!="IgM"].reset_index(drop=True))

    # expected differences in IgM
    outputs.loc[outputs.isotype=="IgM"].reset_index(drop=True).compare(old.loc[old.isotype=="IgM"].reset_index(drop=True))

if __name__=="__main__":
    main()
