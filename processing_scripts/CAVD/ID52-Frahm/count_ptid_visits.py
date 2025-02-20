## -------------------------------------------------------------------------- ##
## Create summaries of Omics data:
## 1. A table with columns [ptid, visit, assay] containing all unique combones
## 2. A second summary table with counts of ptids per tx_group/visit/assay
##
## Author: Beatrix Haddock
## Date: 2025/02/20
## -------------------------------------------------------------------------- ##
import pandas as pd
import numpy as np
import os
## -------------------------------------------------------------------------- ##

# dir to read data in from
parent_dir = '/networks/vtn/lab/SDMC_labscience/studies/VISC/ID52-Frahm/omics_misc/misc_files/data_processing/'

# dir to save outputs to
savedir = '/networks/cavd/Objective 4/GH-VAP/ID52-Frahm/Data/OHSU_omics/'

# read in processed versions of datasets
lcms_metabolomics = pd.read_csv(parent_dir + "/LCMS_Metabolomics/DRAFT_LCMS_Metabolomics_processed_2025-02-19.txt", sep="\t")
gcms_metabolomics = pd.read_csv(parent_dir + "/GCMS_Metabolomics/DRAFT_GCMS_Metabolomics_processed_2025-02-15.txt", sep="\t")
lcms_proteomics = pd.read_csv(parent_dir + "/LCMS_Proteomics/DRAFT_LCMS_Proteomics_processed_2025-02-19.txt", sep="\t")
lcms_lipidomics = pd.read_csv(parent_dir + "/LCMS_Lipidomics/DRAFT_LCMS_Lipidomics_processed_2025-02-20.txt", sep="\t")

# add a sample_type column to proteomics (for subsetting)
lcms_proteomics['sample_type'] = 'Sample'

# concat all data together
datasets = [lcms_metabolomics, gcms_metabolomics, lcms_proteomics, lcms_lipidomics]
all_data = pd.concat([df[['sample_type','ptid','Visit','assay_type']] for df in datasets])

# subset to samples (discard blanks and QCs)
all_data = all_data.loc[all_data.sample_type=="Sample"]
all_data = all_data.drop(columns=['sample_type'])

# subset to unique combinations of ptid/visit/assay
all_data = all_data.drop_duplicates()

# cast ptid (float) to int
all_data.ptid = all_data.ptid.astype(int)

# rename assay column
all_data = all_data.rename(columns={"assay_type":"Assay"})

# save dataset
all_data[['Assay','ptid','Visit']].to_excel(savedir + "assay_ptid_visit_samples.xlsx", index=False)

# read in and merge on treatment assignment
tx_path = '/networks/cavd/Objective 4/GH-VAP/ID52-Frahm/Data/pilot_adata/clinical/TB018_full_rx_pilot.csv'
tx = pd.read_csv(tx_path)
summary = all_data.merge(tx[['ptid','treatment']], on='ptid', how='left')

# groupby assay/tx/visit, and count unique ptids per combo
summary = summary.groupby(['Assay','treatment','Visit'])[['ptid']].count().reset_index()

# cast to wide on visits
summary = pd.pivot_table(summary, index=['Assay','treatment'], columns=['Visit'], values='ptid', fill_value=0)

# save summary dataset
summary.to_excel(savedir + "ohsu_omics_data_summary.xlsx")

# check ptids against clinical dataset
check_clinical_dataset = pd.read_csv("/networks/cavd/obj4/GH-VAP/ID52-Frahm/Data/pilot_adata/clinical/TB018_full_rx_pilot.csv")

# it matches all rows in this dataset where "assay_testing"==0.0
set(check_clinical_dataset.loc[check_clinical_dataset.assay_testing==0.0].ptid).symmetric_difference(all_data.ptid)
