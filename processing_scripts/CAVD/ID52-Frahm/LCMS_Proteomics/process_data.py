## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 03/11/2025
# Purpose: Process LCMS Proteomics data. Create qdata + adata
## ---------------------------------------------------------------------------##

import pandas as pd
import numpy as np
import os
import datetime

def main():
    ## read in data ----------------------------------------------------------##
    datadir = '/networks/cavd/Objective 4/GH-VAP/ID52-Frahm/Data/OHSU_omics/'
    savedir = '/networks/vtn/lab/SDMC_labscience/studies/VISC/ID52-Frahm/omics_misc/misc_files/data_processing/LCMS_Proteomics/'

    proteomics_filepath = datadir + 'GatesM72_PNNL_LCMS_Proteomics_2024_08_FINAL.xlsx'
    peptide_matrices = {}
    for i in range(1,9):
        df = pd.read_excel(proteomics_filepath, sheet_name=f'Peptide_Matrix_P0{i}')
        df['peptide_matrix'] = f'M0{i}'
        peptide_matrices[i] = df

    df_long = {}
    for i in range(1,9):
        df = peptide_matrices[i]
        id_vars = ['Peptide','peptide_matrix'] + [i for i in df.columns if 'empty' in i.lower()] + [i for i in df.columns if 'PooledReference' in i]
        df = df.melt(id_vars=id_vars,
                     value_vars = [i for i in df.columns if 'prepid' in i.lower()],
                     var_name='sample_name',
                     value_name='raw_value'
                    )

        if i < 6:
            df = df.rename(columns = {f'Empty_Plex0{i}': 'Empty_Plex_1',
                                      f'PooledReference_Plex0{i}': 'PooledReference_Plex'})
            df['Empty_Plex_2'] = np.nan
            df['Empty_Plex_3'] = np.nan
        else:
            df = df.rename(columns = {f'Empty_Plex0{i}_1': 'Empty_Plex_1',
                                      f'Empty_Plex0{i}_2': 'Empty_Plex_2',
                                      f'Empty_Plex0{i}_3': 'Empty_Plex_3',
                                      f'PooledReference_Plex0{i}': 'PooledReference_Plex'})

        df_long[i] = df

    df = pd.concat(df_long.values())

    ## merge on sampleid columns----------------------------------------------##
    sampleids = df.sample_name.str.replace("-","_").str.split("_d", expand=True)
    sampleids.columns = ['ptid','Visit']

    # correct labeling issue
    sampleids.loc[sampleids.ptid=='6657_1', 'Visit'] = '37'
    sampleids.loc[sampleids.ptid=='6657_2', 'Visit'] = '0'

    # grab everything from the first column before "_"
    sampleids.ptid = sampleids.ptid.str.split("_", expand=True)[0]

    #grab everything from the second column before "_"
    sampleids['Visit'] = sampleids.Visit.str.split("_", expand=True)[0]

    sampleids['Visit'] = "Day " + sampleids.Visit
    sampleids['sampleid'] = sampleids.ptid + "|" + sampleids.Visit

    df = pd.concat([df, sampleids], axis=1)

    ## merge on metadata -----------------------------------------------------##
    peptide_metadata = pd.read_excel(proteomics_filepath, sheet_name = 0)
    df = df.merge(peptide_metadata, on='Peptide', how='left')

    timestamp = datetime.datetime.fromtimestamp(os.path.getmtime(proteomics_filepath)).date().isoformat()
    sdmc_processing_datetime = datetime.datetime.now().replace(microsecond=0).isoformat()

    std_metadata = pd.DataFrame({
        'network': ['GHDC'],
        'protocol': ['ID52-Frahm'],
        'assay_lab_name': ['Tafesse OHSU/PNNL'],
        'assay_type': ['GCMS Proteomics'],
        # 'transform': ['log'],
        'Unit': ['intensity'],
        'upload_date': [timestamp],
        'filename': [proteomics_filepath.rpartition("/")[-1]],
        # 'LLOD': [0.0],
        'sdmc_processing_datetime':[sdmc_processing_datetime],
    })

    df = df.merge(std_metadata, how='cross')

    ## reorder ---------------------------------------------------------------##
    reorder = [
        'network',
        'protocol',
        'assay_lab_name',
        'assay_type',
        'ptid',
        'Visit',
        'sampleid',
        'sample_name',
        'Peptide',
        'peptide_matrix',
        # 'bioid',
        'raw_value',
        # 'Value',
        'PooledReference_Plex',
        'Empty_Plex_1',
        'Empty_Plex_2',
        'Empty_Plex_3',
        # 'LLOD',
        'Protein',
        'Spectra',
        'MSGF_SpecProb',
        # 'transform',
        'Unit',
        'upload_date',
        'filename',
        'sdmc_processing_datetime',
    ]

    check = set(df.columns).symmetric_difference(reorder)
    if len(check) > 0:
        raise Warning("trying to add or drop columns")
    df = df[reorder]


    ## save to txt -----------------------------------------------------------##
    today = datetime.date.today().isoformat()
    df.to_csv(savedir + f"DRAFT_LCMS_Proteomics_processed_{today}.txt", sep="\t", index=False)

    ## pivot summary ---------------------------------------------------------##
    # summary = pd.pivot_table(df,
    #                index='Peptide',
    #                columns=['ptid','Visit'],
    #                aggfunc='count',
    #                fill_value=0
    #               )[['raw_value']]
    # summary.to_csv(savedir + "LCMS_Proteomics_pivot_summary.txt", sep="\t")

    ## create adata ----------------------------------------------------------##

    # make sure that these two columns create a unique row identifier
    df[['Peptide','sampleid']].drop_duplicates().shape[0] == len(df)

    # create a dataset with unique row for each peptide/sampleid column.
    # each of these is associated with a bunch of metadata
    # note 'peptide_matrix' refers to the tab of the input data workbook each row of data comes from, so will be 'NA' for all new combos of sampleid/peptide
    adata = df[['Peptide','Protein', 'Spectra', 'MSGF_SpecProb',]].drop_duplicates().merge(
        df[['network',
            'protocol',
            'assay_lab_name',
            'assay_type',
            'ptid',
            'Visit',
            'sampleid',
            'sample_name',
            'Unit',
            'upload_date',
            'filename',
            'sdmc_processing_datetime']].drop_duplicates(),
        how='cross'
    )

    # add bioid column. peptide is unique here.
    adata['bioid'] = adata['Peptide']

    # merge data onto this frame of peptide/sampleid pairs. keep NA for all rows where we didnt get anything from the lab
    adata = adata.merge(df[['Peptide','sampleid','peptide_matrix','raw_value','PooledReference_Plex']], on=['Peptide','sampleid'], how='left')

    # rename for cross-dataset consistency
    adata = adata.rename(columns={
        'Peptide': 'Antigen',
    })

    # calculate primary 'value' column following jennifer kyle's instructions
    adata['Value'] = adata.raw_value / adata.PooledReference_Plex

    # add metadata for programming
    adata_metadata = pd.DataFrame({
        'LLOD': [0.0],
        'transform': ['log'],
    })
    adata = adata.merge(adata_metadata, how='cross')

    # reorder columns
    reorder_adata = [
        # 'network',
        # 'protocol',
        # 'assay_lab_name',
        # 'assay_type',
        'ptid',
        'Visit',
        'sampleid',
        'sample_name',
        'Antigen',
        'peptide_matrix',
        'bioid',
        'raw_value',
        'Value',
        'PooledReference_Plex',
        # 'Empty_Plex_1',
        # 'Empty_Plex_2',
        # 'Empty_Plex_3',
        'LLOD',
        'Protein',
        'Spectra',
        'MSGF_SpecProb',
        'transform',
        'Unit',
        'upload_date',
        'filename',
        'sdmc_processing_datetime',
    ]

    check = set(adata.columns).symmetric_difference(reorder_adata)
    if len(check)>0:
        print(f"adding or dropping columns for adata: {check}")

    adata = adata[reorder_adata]

    # save to .txt file
    today = datetime.date.today().isoformat()
    adata.to_csv(savedir + f'ohsu_lcms_proteomics_adata_{today}.txt', sep="\t", index=False)

if __name__=="__main__":
    main()
