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
    sampleids = lcms_lipid.sample_name.str[len("Gates_Mtb_"):].str.split("_d", expand=True)

    sampleids.columns = ['ptid','Visit']

    sampleids.loc[lcms_lipid.sample_type!='Sample','ptid'] = np.nan
    sampleids.loc[lcms_lipid.sample_type!='Sample','Visit'] = np.nan

    sampleids['Visit'] = sampleids.Visit.str.split("_", expand=True)[0]

    sampleids.Visit.unique()

    # correct labeling issue
    sampleids.loc[sampleids.ptid=='6657','Visit'] = sampleids.loc[sampleids.ptid=='6657'].Visit.map({'0-1':'37','0-2':'0'})
    sampleids['Visit'] = "Day " + sampleids.Visit

    sampleids['sampleid'] = sampleids['ptid'] + "|" + sampleids['Visit']

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

    adata['bioid'] = adata[['Antigen','RT','Adduct type']].astype(str).apply('|'.join, axis=1)
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
    summary = pd.pivot_table(adata.loc[adata.sampleid.notna()],
                             index=['Antigen','RT','Adduct type'],
                              columns='Visit',
                             aggfunc='count',
                             fill_value=0)[['Value']]
    summary = summary.rename(columns={'Value':'ptid_count'})
    summary.to_excel(savedir + "LCMS_Lipidomics_bioid_Visit_summary.xlsx")

    # Antigen','RT','Adduct type','sample_name
    summary2 = pd.pivot_table(adata.loc[adata.sampleid.notna()],
                              index=['Antigen','RT','Adduct type'],
                              columns=['ptid','Visit'],
                              aggfunc='count',
                              fill_value=0)[['Value']]
    summary2 = summary2.rename(columns={'Value':'ptid_count'})
    summary2.to_excel(savedir + "LCMS_Lipidomics_bioid_ptid_Visit_summary.xlsx")

    ## checks ----------------------------------------------------------------##
    check_clinical_dataset = pd.read_csv("/networks/cavd/obj4/GH-VAP/ID52-Frahm/Data/pilot_adata/clinical/TB018_full_rx_pilot.csv")
    set(lcms_lipid.ptid.dropna().astype(int)).difference(check_clinical_dataset.ptid.astype(int))

if __name__=="__main__":
    main()
