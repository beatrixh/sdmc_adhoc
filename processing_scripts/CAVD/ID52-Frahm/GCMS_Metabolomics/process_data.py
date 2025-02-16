## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 02/15/2025
# Purpose: Process GCMS Metabolomics data. Create qdata + adata
## ---------------------------------------------------------------------------##
import pandas as pd
import numpy as np
import os
import datetime

def main():
    datadir = '/networks/cavd/Objective 4/GH-VAP/ID52-Frahm/Data/OHSU_omics/'
    savedir = '/networks/vtn/lab/SDMC_labscience/studies/VISC/ID52-Frahm/omics_misc/misc_files/data_processing/GCMS_Metabolomics/'

    # read in first sheet (primary data)
    gcms_metab_filepath = datadir + "GatesM72_PNNL_GCMS_Metabolomics_2024_08_FINAL.xlsx"
    gcms_metab = pd.read_excel(gcms_metab_filepath, sheet_name = 0, header=None)

    # combine top two rows into one piece of info (one per sample?) as column names
    gcms_metab.columns = [f"{j};{i}" for (i,j) in zip(gcms_metab.iloc[0],gcms_metab.iloc[1])]
    gcms_metab = gcms_metab.iloc[2:]

    # these are the same; so Metabolite a unique identifier
    gcms_metab['Metabolite;nan'].nunique(), len(gcms_metab)

    # cast to long
    gcms_metab = gcms_metab.iloc[2:].melt(id_vars=['Metabolite;nan','Annotation description;nan'],
                                            value_vars=list(gcms_metab.columns[2:]),
                                            value_name='Value'
                                           )

    # separate top two rows back into two pieces of info (two columns now)
    gcms_metab = pd.concat([gcms_metab, gcms_metab.variable.str.split(";", expand=True)], axis=1)
    gcms_metab = gcms_metab.drop(columns='variable')

    # filter out rows with weird (non numeric) values
    gcms_metab['qc_filter'] = False
    gcms_metab.loc[gcms_metab.Value=='tegration parameters!', 'qc_filter'] = True
    gcms_metab.loc[gcms_metab.Value=='tegration parameters!', 'Value'] = np.nan

    # rename columns
    gcms_metab.columns = ['Antigen', 'metab_annotation', 'Value', 'sample_name', 'sample_type', 'qc_filter']
    gcms_metab[['Antigen','sample_name']].drop_duplicates().shape, gcms_metab.shape

    # parse out ptid, Visit, sampleid, then merge on
    sampleids = gcms_metab.sample_name.str[len("Gates_Mtb_"):].str.replace("-","_").str.split("_", expand=True).iloc[:,:2]
    sampleids.columns = ['ptid','Visit']
    sampleids['Visit'] = "Day " + sampleids.Visit.str[1:]
    sampleids.loc[gcms_metab.sample_type!='Sample','ptid'] = np.nan
    sampleids.loc[gcms_metab.sample_type!='Sample','Visit'] = np.nan

    sampleids['sampleid'] = sampleids['ptid'] + "|" + sampleids['Visit']
    gcms_metab = pd.concat([gcms_metab, sampleids], axis=1)

    # verify these create a bioid
    bioid_is_unique = len(gcms_metab[['Antigen','sample_name']].drop_duplicates()) == len(gcms_metab)
    if not bioid_is_unique:
        raise Exception("Bioid not unique")

    # merge on bioid
    gcms_metab['bioid'] = gcms_metab['Antigen'] + "|" + gcms_metab['sample_name']

    timestamp = datetime.datetime.fromtimestamp(os.path.getmtime(gcms_metab_filepath)).date().isoformat()

    gcms_metab_metadata = pd.DataFrame({
        'network': ['GHDC'],
        'protocol': ['ID52-Frahm'],
        'lab_name': ['Tafesse OHSU/PNNL'],
        'assay_type': ['GCMS Metabolomics'],
        'transform': ['log'],
        'Unit': ['Relative peak area'],
        'upload_date': [timestamp],
        'filename': [gcms_metab_filepath.rpartition("/")[-1]],
        'LLOD': [0.0]
    })

    gcms_metab = gcms_metab.merge(gcms_metab_metadata, how='cross')

    reorder = [
        'network',
        'protocol',
        'lab_name',
        'assay_type',
        'ptid',
        'Visit',
        'sampleid',
        'Antigen',
        'bioid',
        'Value',
        'LLOD',
        'sample_name',
        'sample_type',
        'transform',
        'Unit',
        'upload_date',
        'filename',
        'metab_annotation',
        'qc_filter'
    ]

    # ensure this doesn't try to add or drop columns
    mismatch = set(gcms_metab.columns).symmetric_difference(reorder)
    if len(mismatch) > 0:
        raise Warning(f"Trying to add or drop colummns: {mismatch}")

    gcms_metab = gcms_metab[reorder]

    ## Checks ----------------------------------------------------------------##

    # every sample has 135 antigens
    gcms_metab.groupby(['sample_name']).Antigen.nunique().value_counts()

    # every antigen has 164 samples
    gcms_metab.groupby(['Antigen']).sample_name.nunique().value_counts()

    # clinical dataset from andrew
    check_clinical_dataset = pd.read_csv("/networks/cavd/obj4/GH-VAP/ID52-Frahm/Data/pilot_adata/clinical/TB018_full_rx_pilot.csv")

    # exact match between data ptids and assay_testing==0.0 rows
    set(gcms_metab.ptid.dropna().astype(int)).symmetric_difference(check_clinical_dataset.loc[check_clinical_dataset.assay_testing==0.0].ptid.astype(int))

    # check against other omics datasets
    set(gcms_metab.sampleid.dropna()).difference(lcms_metab.sampleid.dropna())
    set(lcms_metab.sampleid.dropna()).difference(gcms_metab.sampleid.dropna())

    # missing 6657|Day 37
    set(gcms_metab.sampleid.dropna()).difference(prot.sampleid.dropna())
    set(prot.sampleid.dropna()).difference(gcms_metab.sampleid.dropna())

    set(gcms_metab.sampleid.dropna()).difference(lipid.sampleid.dropna())
    set(lipid.sampleid.dropna()).difference(gcms_metab.sampleid.dropna())

    ## save to txt -----------------------------------------------------------##
    today = datetime.date.today().isoformat()
    gcms_metab.to_csv(savedir + f"DRAFT_GCMS_Metabolomics_processed_{today}.txt", sep="\t", index=False)

    ## pivot summary ---------------------------------------------------------##
    pivot_summary = gcms_metab.copy()
    pivot_summary[['ptid','Visit','Antigen']] = pivot_summary[['ptid','Visit','Antigen']].fillna("NA")
    summary = pd.pivot_table(pivot_summary,
                             index=['ptid','Visit'],
                             columns=['Antigen'],
                             aggfunc='count',
                             fill_value=0
                            )[['Value']]
    summary.to_csv(savedir + "GCMS_Metabolomics_pivot_summary.txt", sep="\t")

if __name__=="__main__":
    main()
