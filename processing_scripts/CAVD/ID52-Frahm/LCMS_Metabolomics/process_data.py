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
    ## -----------------------------------------------------------------------##
    datadir = '/networks/cavd/Objective 4/GH-VAP/ID52-Frahm/Data/OHSU_omics/'
    savedir = '/networks/vtn/lab/SDMC_labscience/studies/VISC/ID52-Frahm/omics_misc/misc_files/data_processing/LCMS_Metabolomics/'

    ## read in data ----------------------------------------------------------##
    lcms_metab_filepath = datadir + 'GatesM72_PNNL_LCMS_Metabolomics_2024_08_FINAL.xlsx'

    lcms_metab = {}
    tabs = ['Tab 2 - RP Pos', 'Tab 3 - RP Neg', 'Tab 4 - HILIC Pos', 'Tab 5 - HILIC Neg']
    for pm in tabs:
        s = pd.read_excel(datadir + 'GatesM72_PNNL_LCMS_Metabolomics_2024_08_FINAL.xlsx', sheet_name=pm, header=None)

        # move first two rows into column names
        s.columns = [f"{j};{i}" if pd.notna(i) else f"{j}" for (i,j) in zip(s.iloc[0],s.iloc[1])]
        s = s.iloc[2:].reset_index(drop=True)

        # cast to long
        s = s.melt(
            id_vars=list(s.columns)[:12],
            value_vars=list(s.columns)[12:],
            value_name='Value'
        )

        # split formerly first two rows into two separate columns
        s = pd.concat([
            s.drop(columns='variable'),
            s.variable.str.split(";", expand=True).rename(columns={0:'sample_name', 1:'sample_type'})
        ], axis=1)
        phase_mode = pm.split(" - ")[-1]
        s['phase_mode'] = phase_mode
        lcms_metab[phase_mode] = s

    # confirmed they all have the same columns
    for i in lcms_metab.keys():
        for j in lcms_metab.keys():
            if i!=j:
                diff = set(lcms_metab[i].columns).symmetric_difference(lcms_metab[j].columns)
                if len(diff) > 0:
                    print(i, j)

    lcms_metab = pd.concat(lcms_metab.values())

    ## merge on sample id columns --------------------------------------------##
    sampleids = lcms_metab['sample_name'].str[len('Gates_Mtb_'):].str.split("_", expand=True)
    sampleids = sampleids.iloc[:,:2]

    sampleids.columns = ['ptid','Visit']

    sampleids.loc[lcms_metab.sample_type!='Sample', 'ptid'] = np.nan
    sampleids.loc[lcms_metab.sample_type!='Sample', 'Visit'] = np.nan

    sampleids.Visit = sampleids.Visit.str.split("-", expand=True)[0]

    sampleids['Visit'] = "Day " + sampleids['Visit'].str[1:]
    sampleids['sampleid'] = sampleids.ptid + "|" + sampleids.Visit

    lcms_metab = pd.concat([lcms_metab, sampleids], axis=1)

    ## merge on metadata -----------------------------------------------------##
    timestamp = datetime.datetime.fromtimestamp(os.path.getmtime(lcms_metab_filepath)).date().isoformat()
    sdmc_processing_datetime = datetime.datetime.now().replace(microsecond=0).isoformat()

    lcms_metab_metadata = pd.DataFrame({
        'network': ['GHDC'],
        'protocol': ['ID52-Frahm'],
        'assay_lab_name': ['Tafesse OHSU/PNNL'],
        'assay_type': ['LCMS Metabolomics'],
        'Unit': ['Relative peak area'],
        'upload_date': [timestamp],
        'filename': [lcms_metab_filepath.rpartition("/")[-1]],
        'sdmc_processing_datetime':[sdmc_processing_datetime]
    })
    lcms_metab = lcms_metab.merge(lcms_metab_metadata, how='cross')

    # trim leading/trailing whitespace
    lcms_metab.columns = [i.strip() for i in lcms_metab.columns]

    ## reorder columns -------------------------------------------------------##
    # I typically put the sample columns first (ptid, Visit, sampleid), then the assay columns (eg Antigen, Isotype),
    # then Value and other metadata columns like units and transform and additional columns.

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
        'phase_mode',
        'Value',
        'm/z',
        'RT [min]',
        'Standardized name',
        'Super class',
        'Main class',
        'Sub class',
        'Formula',
        'Annotation MW',
        'Annot. DeltaMass [ppm]',
        'Reference Ion',
        'Unit',
        'upload_date',
        'filename',
        'sdmc_processing_datetime',
    ]

    set(lcms_metab.columns).symmetric_difference(reorder_qdata)

    lcms_metab = lcms_metab[reorder_qdata]

    ## save qdata ------------------------------------------------------------##
    today = datetime.date.today().isoformat()
    lcms_metab.to_csv(savedir + f"DRAFT_LCMS_Metabolomics_processed_{today}.txt", sep="\t", index=False)

    ## adata -----------------------------------------------------------------##

    adata = lcms_metab.copy()
    adata = adata.rename(columns={
        'Metabolite': 'Antigen',
        'Annotation description': 'metab_annotation',
        'RT [min]': 'RT'
    })

    adata_metadata = pd.DataFrame({
        'LLOD': [0.0],
        'transform': ['log'],
    })

    adata = adata.merge(adata_metadata, how='cross')

    adata['bioid'] = adata[['Antigen','RT', 'phase_mode']].astype(str).apply('|'.join, axis=1)
    len(adata[['bioid','sample_name']].drop_duplicates()) == len(adata)

    reorder_adata = [
        'ptid',
        'Visit',
        'sampleid',
        'sample_name',
        'sample_type',
        'Antigen',
        'metab_annotation',
        'bioid',
        'phase_mode',
        'Value',
        'LLOD',
        'm/z',
        'RT',
        'Standardized name',
        'Super class',
        'Main class',
        'Sub class',
        'Formula',
        'Annotation MW',
        'Annot. DeltaMass [ppm]',
        'Reference Ion',
        'transform',
        'Unit',
        'upload_date',
        'filename',
    ]

    set(reorder_adata).symmetric_difference(adata.columns)

    adata = adata[reorder_adata]

    ## save adata ------------------------------------------------------------##

    # I think you were going to save something like qdata first and then produce adata right? The qdata can go in the same folder as the raw data: /networks/cavd/obj4/GH-VAP/ID52-Frahm/Data/OHSU_omics whereas the adata sets can be placed here: networks/cavd/obj4/GH-VAP/ID52-Frahm/Data/pilot_adata (use the naming convention of other files in the folder with a prefix and date. Use prefix ohsu_XXXX for whichever assay it is.
    adata.to_csv(savedir + f'ohsu_lcms_metabolomics_adata_{today}.txt', sep="\t", index=False)

    ## checks ----------------------------------------------------------------##
    check_clinical_dataset = pd.read_csv("/networks/cavd/obj4/GH-VAP/ID52-Frahm/Data/pilot_adata/clinical/TB018_full_rx_pilot.csv")
    set(lcms_metab.ptid.dropna().astype(int)).difference(check_clinical_dataset.ptid.astype(int))

    ## pivot summaries -------------------------------------------------------##
    summary = adata.copy()
    summary[['Antigen','RT', 'phase_mode','sample_name']] = summary[['Antigen','RT', 'phase_mode','sample_name']].fillna("NA")

    summary_with_RT = pd.pivot_table(summary,
                   index=['Antigen','RT','phase_mode'],
                   columns='sample_name',
                   aggfunc='count',
                   fill_value=0
                  )[['Value']]
    summary_with_RT.to_csv(savedir + "LCMS_Metabolomics_pivot_summary_with_RT.txt", sep="\t")

    summary_without_RT = pd.pivot_table(summary,
                   index=['Antigen','phase_mode'],
                   columns='sample_name',
                   aggfunc='count',
                   fill_value=0
                  )[['Value']]
    summary_without_RT.to_csv(savedir + "LCMS_Metabolomics_pivot_summary_with_Antigen_only.txt", sep="\t")


if __name__=="__main__":
    main()
