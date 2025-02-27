## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 02/25/2025
# Purpose: Process LCMS Metabolomics data. Create qdata + adata
## ---------------------------------------------------------------------------##
import pandas as pd
import numpy as np
import os
import datetime

def main():
    ## read in data ----------------------------------------------------------##
    datadir = '/networks/cavd/Objective 4/GH-VAP/ID52-Frahm/Data/OHSU_omics/'
    savedir = '/networks/vtn/lab/SDMC_labscience/studies/VISC/ID52-Frahm/omics_misc/misc_files/data_processing/LCMS_Lipidomics/'

    lipid_filename = datadir + 'GatesM72_PNNL_LCMS_Lipidomics_2024_08_FINAL.xlsx'
    lcms_lipid = {}
    for tab in ['Tab 1 - Pos', 'Tab 2- Neg']:
        df = pd.read_excel(datadir + 'GatesM72_PNNL_LCMS_Lipidomics_2024_08_FINAL.xlsx', sheet_name=tab, header=None)

        # move first two rows into column names
        df.columns = [f"{j};{i}" if pd.notna(i) else f"{j}" for (i,j) in zip(df.iloc[0],df.iloc[1])]
        df = df.iloc[2:].reset_index(drop=True)

        #cast to long
        df = df.melt(
            id_vars=list(df.columns)[:10],
            value_vars=list(df.columns)[10:],
            value_name='Value'
        )

        # split formerly first two rows into two separate columns
        df = pd.concat([
            df.drop(columns='variable'),
            df.variable.str.split(";", expand=True).rename(columns={0:'sample_name', 1:'sample_type'})
        ], axis=1)

        # adding info following 'M72_ID52_OHSU_assay_questions_19Aug2024-Kyle.docx'
        df['phase_mode'] = "RP " + tab.split("-")[-1].strip()

        lcms_lipid[tab] = df

    lcms_lipid = pd.concat(lcms_lipid.values()).reset_index(drop=True)

    ## merge on sample ids ---------------------------------------------------##
    sampleids = lcms_lipid.loc[lcms_lipid.sample_type=="Sample"].sample_name.str[len("Gates_Mtb_"):].str.split("_", expand=True).iloc[:,:2]
    sampleids[1] = sampleids[1].str.split("-", expand=True)[0]

    sampleids.columns = ['ptid', 'Visit']

    sampleids['Visit'] = "Day " + sampleids.Visit.str[1:]
    sampleids['sampleid'] = sampleids.ptid + "|" + sampleids.Visit

    lcms_lipid = pd.concat([lcms_lipid, sampleids], axis=1)

    ## merge on metadata -----------------------------------------------------##\
    timestamp = datetime.datetime.fromtimestamp(os.path.getmtime(lipid_filename)).date().isoformat()
    sdmc_processing_datetime = datetime.datetime.now().replace(microsecond=0).isoformat()

    lipid_metadata = pd.DataFrame({
        'network': ['GHDC'],
        'protocol': ['ID52-Frahm'],
        'assay_lab_name': ['Tafesse OHSU/PNNL'],
        'assay_type': ['LCMS Lipidomics'],
        # 'transform': ['log'],
        'Unit': ['Relative peak apex intensity'],
        'upload_date': [timestamp],
        'filename': [lipid_filename.rpartition("/")[-1]],
        'sdmc_processing_datetime':[sdmc_processing_datetime],
        # 'LLOD': [0.0]
    })

    lcms_lipid = lcms_lipid.merge(lipid_metadata, how = 'cross')

    ## reorder data ----------------------------------------------------------##
    reorder_qdata = [
        'network',
        'protocol',
        'assay_lab_name',
        'assay_type',
        'ptid',
        'Visit',
        'sampleid',
        'sample_name',
        'sample_type',
        'Lipid',
        'phase_mode',
        # 'bioid',
        'Value',
        'S/N average',
        'MS/MS spectrum',
        'Average Mz',
        'Average Rt(min)',
        'Adduct type',
        'Reference m/z',
        'Annot. DeltaMass [ppm]',
        'Formula',
        'Ontology',
        # 'transform',
        # 'LLOD',
        'Unit',
        'upload_date',
        'filename',
        'sdmc_processing_datetime',
    ]

    set(lcms_lipid.columns).symmetric_difference(reorder_qdata)

    lcms_lipid = lcms_lipid[reorder_qdata]

    ## save to txt -----------------------------------------------------------##

    today = datetime.date.today().isoformat()
    lcms_lipid.to_csv(savedir + f"DRAFT_LCMS_Lipidomics_processed_{today}.txt", sep="\t", index=False)

    ## adata -----------------------------------------------------------------##
    adata = lcms_lipid.copy()
    adata = adata.rename(columns={
        'Lipid': 'Antigen',
        'Average Rt(min)': 'RT'
    })

    adata_metadata = pd.DataFrame({
        'LLOD': [0.0],
        'transform': ['log'],
    })

    adata = adata.merge(adata_metadata, how='cross')

    adata['bioid'] = adata[['Antigen','RT','Adduct type','sample_name']].astype(str).apply('|'.join, axis=1)
    adata[['bioid','sample_name']].drop_duplicates().shape[0] == len(adata)

    reorder_adata = [
        'ptid',
        'Visit',
        'sampleid',
        'sample_name',
        'sample_type',
        'Antigen',
        'phase_mode',
        'bioid',
        'Value',
        'S/N average',
        'MS/MS spectrum',
        'Average Mz',
        'RT',
        'Adduct type',
        'Reference m/z',
        'Annot. DeltaMass [ppm]',
        'Formula',
        'Ontology',
        'transform',
        'LLOD',
        'Unit',
        'upload_date',
        'filename',
    ]

    set(adata.columns).symmetric_difference(reorder_adata)
    adata = adata[reorder_adata]

    adata.to_csv(savedir + f'ohsu_lcms_lipidolomics_adata_{today}.txt', sep="\t", index=False)

    ## pivot summary ---------------------------------------------------------##
    summary = adata.copy()
    summary[['Antigen','S/N average','sample_name','phase_mode']] = summary[['Antigen','S/N average','sample_name','phase_mode']].fillna("NA")

    pivot_summary = pd.pivot_table(
        summary,
        index=['Antigen','S/N average','phase_mode'],
        columns=['ptid','Visit'],
        aggfunc='count',
        fill_value=0
    )[['Value']]

    pivot_summary.to_csv(savedir + "LCMS_Lipidomics_pivot_summary_with_sn_average.txt", sep="\t")

    pivot_summary_just_antigen = pd.pivot_table(
        summary,
        index=['Antigen','phase_mode'],
        columns=['ptid','Visit'],
        aggfunc='count',
        fill_value=0
    )[['Value']]

    pivot_summary_just_antigen.to_csv(savedir + "LCMS_Lipidomics_pivot_summary_without_sn_average.txt", sep="\t")

    ## checks ----------------------------------------------------------------##
    check_clinical_dataset = pd.read_csv("/networks/cavd/obj4/GH-VAP/ID52-Frahm/Data/pilot_adata/clinical/TB018_full_rx_pilot.csv")
    set(lcms_lipid.ptid.dropna().astype(int)).difference(check_clinical_dataset.ptid.astype(int))

if __name__=="__main__":
    main()
