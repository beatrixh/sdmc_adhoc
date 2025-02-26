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
    ## -----------------------------------------------------------------------##
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

    gcms_metab['qc_filter'] = False
    gcms_metab.loc[gcms_metab.Value=='tegration parameters!', 'qc_filter'] = True
    gcms_metab.loc[gcms_metab.Value=='tegration parameters!', 'Value'] = np.nan

    gcms_metab = gcms_metab.rename(columns={
        'Metabolite;nan':'Metabolite',
        'Annotation description;nan':'Annotation description',
        0:'sample_name',
        1:'sample_type'
    })

    # merge on sample id columns
    sampleids = gcms_metab.sample_name.str[len("Gates_Mtb_"):].str.replace("-","_").str.split("_", expand=True).iloc[:,:2]
    sampleids.columns = ['ptid','Visit']
    sampleids['Visit'] = "Day " + sampleids.Visit.str[1:]

    sampleids.loc[gcms_metab.sample_type!='Sample','ptid'] = np.nan
    sampleids.loc[gcms_metab.sample_type!='Sample','Visit'] = np.nan

    sampleids['sampleid'] = sampleids['ptid'] + "|" + sampleids['Visit']
    gcms_metab = pd.concat([gcms_metab, sampleids], axis=1)

    timestamp = datetime.datetime.fromtimestamp(os.path.getmtime(gcms_metab_filepath)).date().isoformat()
    sdmc_processing_datetime = datetime.datetime.now().replace(microsecond=0).isoformat()

    gcms_metab_metadata = pd.DataFrame({
        'network': ['GHDC'],
        'protocol': ['ID52-Frahm'],
        'assay_lab_name': ['Tafesse OHSU/PNNL'],
        'assay_type': ['GCMS Metabolomics'],
        'Unit': ['Relative peak area'],
        'upload_date': [timestamp],
        'filename': [gcms_metab_filepath.rpartition("/")[-1]],
        'sdmc_processing_datetime':[sdmc_processing_datetime]
    })

    gcms_metab = gcms_metab.merge(gcms_metab_metadata, how='cross')

    # From andrew: I typically put the sample columns first (ptid, Visit, sampleid), then the assay columns (eg Antigen, Isotype), then Value and other metadata columns like units and transform and additional columns.

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
        'Metabolite',
        'Annotation description',
        'Value',
        'Unit',
        'qc_filter',
        'upload_date',
        'filename',
        'sdmc_processing_datetime'
    ]

    # check not dropping or adding columns
    check = set(reorder_qdata).symmetric_difference(gcms_metab.columns)
    if len(check)>0:
        raise Warning(f"trying to add or drop columns {check}")
    gcms_metab = gcms_metab[reorder_qdata]

    ## save qdata to txt -----------------------------------------------------##
    today = datetime.date.today().isoformat()
    gcms_metab.to_csv(savedir + f"DRAFT_GCMS_Metabolomics_processed_{today}.txt", sep="\t", index=False)

    ## adata -----------------------------------------------------------------##
    adata = gcms_metab.copy()
    adata = adata.rename(columns={
        'Metabolite': 'Antigen',
        'Annotation description': 'metab_annotation',
    })

    adata_metadata = pd.DataFrame({
        'LLOD': [0.0],
        'transform': ['log'],
    })
    adata = adata.merge(adata_metadata, how='cross')

    adata['bioid'] = adata['Antigen']

    check = len(adata[['bioid','sample_name']].drop_duplicates()) == len(adata)
    if not check:
        raise Warning("Bioid not unique")

    reorder_adata = [
        'ptid',
        'Visit',
        'sampleid',
        'sample_name',
        'sample_type',
        'Antigen',
        'metab_annotation',
        'bioid',
        'Value',
        'Unit',
        'LLOD',
        'transform',
        'qc_filter',
        'upload_date',
        'filename',
    ]

    check = set(reorder_adata).symmetric_difference(adata.columns)
    if len(check)>0:
        raise Warning(f"trying to add or drop columns {check}")
    adata = adata[reorder_adata]

    ## save adata to txt -----------------------------------------------------##
    # I think you were going to save something like qdata first and then produce adata right? The qdata can go in the same folder as the raw data: /networks/cavd/obj4/GH-VAP/ID52-Frahm/Data/OHSU_omics whereas the adata sets can be placed here: networks/cavd/obj4/GH-VAP/ID52-Frahm/Data/pilot_adata (use the naming convention of other files in the folder with a prefix and date. Use prefix ohsu_XXXX for whichever assay it is.
    adata.to_csv(savedir + f'ohsu_gcms_metabolomics_adata_{today}', sep="\t", index=False)

    ## Checks ----------------------------------------------------------------##

    # every sample has 135 antigens
    adata.groupby(['sample_name']).Antigen.nunique().value_counts()

    # every antigen has 164 samples
    adata.groupby(['Antigen']).sample_name.nunique().value_counts()

    # clinical dataset from andrew
    check_clinical_dataset = pd.read_csv("/networks/cavd/obj4/GH-VAP/ID52-Frahm/Data/pilot_adata/clinical/TB018_full_rx_pilot.csv")

    # exact match between data ptids and assay_testing==0.0 rows
    set(adata.ptid.dropna().astype(int)).symmetric_difference(check_clinical_dataset.loc[check_clinical_dataset.assay_testing==0.0].ptid.astype(int))

    # check against other omics datasets
    path = '/networks/vtn/lab/SDMC_labscience/studies/VISC/ID52-Frahm/omics_misc/misc_files/data_processing/LCMS_Metabolomics/DRAFT_LCMS_Metabolomics_processed_2025-02-15.txt'
    lcms_metab = pd.read_csv(path, sep="\t")

    path = '/networks/vtn/lab/SDMC_labscience/studies/VISC/ID52-Frahm/omics_misc/misc_files/data_processing/LCMS_Proteomics/DRAFT_LCMS_Proteomics_processed_2025-02-15.txt'
    prot = pd.read_csv(path, sep="\t")

    path = '/networks/vtn/lab/SDMC_labscience/studies/VISC/ID52-Frahm/omics_misc/misc_files/data_processing/LCMS_Lipidomics/DRAFT_LCMS_Lipidomics_processed_2025-02-15.txt'
    lipid = pd.read_csv(path, sep="\t")

    set(adata.sampleid.dropna()).difference(lcms_metab.sampleid.dropna())
    set(lcms_metab.sampleid.dropna()).difference(adata.sampleid.dropna())

    # missing 6657|Day 37
    set(adata.sampleid.dropna()).difference(prot.sampleid.dropna())
    set(prot.sampleid.dropna()).difference(adata.sampleid.dropna())

    set(adata.sampleid.dropna()).difference(lipid.sampleid.dropna())
    set(lipid.sampleid.dropna()).difference(adata.sampleid.dropna())

    ## pivot summary ---------------------------------------------------------##
    pivot_summary = adata.copy()
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
