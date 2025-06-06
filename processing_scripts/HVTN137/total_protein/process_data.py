## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 06/05/2025
# Purpose:  - Process HVTN137 total protein data from tomaras lab
#           - Generate pivot summary
## ---------------------------------------------------------------------------##
import pandas as pd
import numpy as np
import datetime
import yaml
import os

import sdmc_tools.process as sdmc
import sdmc_tools.constants as constants
import sdmc_tools.access_ldms as access_ldms

def main():
    ## standard processing \\\\-----------------------------------------------##
    input_data_path = "/trials/vaccine/p137/s001/qdata/LabData/BAMA_Mucosal/total_protein_pass-through/HVTN137_Mucosal_Total_Protein_20250331_A.txt"
    data = pd.read_csv(input_data_path, sep="\t")

    data = data.rename(columns={
        "global_identifier":"guspec",
        "observed_concentration":"concentration"
        })

    ldms = access_ldms.pull_one_protocol('HVTN', 137)
    ldms = ldms.loc[ldms.guspec.isin(data.guspec)].drop_duplicates()

    for col in ['drawdm', 'drawdd', 'drawdy', 'lstudy']:
        ldms[col] = ldms[col].astype(str).astype(float).astype(int)

    metadata_dict = {
        "network": "HVTN",
        "upload_lab_id": "GT",
        "assay_lab_name": "Tomaras Lab (Duke)",
        "assay_kit": "Quant-iT",
        "assay_type": "Total Protein",
        "specrole": "Sample",
    }

    outputs = sdmc.standard_processing(
        input_data=data,
        input_data_path=input_data_path,
        guspec_col="guspec",
        network='hvtn',
        metadata_dict=metadata_dict,
        ldms=ldms
    )

    if outputs.spectype.isna().sum() > 0:
        raise Exception("THERE ARE NANS IN SPEDCTYPE; NEED TO FIX")

    ## check these are equivalent before dropping ----------------------------##
    checks = 0
    checks += (outputs.visit_identifier.astype(float) != outputs.visitno.astype(float)).sum()
    checks += (outputs.sample_identifier.astype(float) != outputs.ptid.astype(float)).sum()
    assert checks == 0

    outputs = outputs.drop(columns=[
        "assay_platform", "lab", "study", "visit_identifier",
        "timepoint", "timepoint_units", "isotype", "sample_identifier"
    ])

    ## save to txt -----------------------------------------------------------##
    savedir = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN137/assays/total_protein/misc_files/data_processing/'
    today = datetime.date.today().isoformat()
    savepath = savedir + f'HVTN137_total_protein_processed_{today}.txt'
    outputs.to_csv(savepath, index=False, sep="\t")

    ## pivot summary ---------------------------------------------------------##
    pivot_summary = pd.pivot_table(outputs[["ptid","visitno","spectype","concentration"]],
                                   index=["ptid", "visitno"],
                                   columns="spectype",
                                   aggfunc='count',
                                   dropna=False,
                                   fill_value=0
                                  )

    pivot_summary = pivot_summary.droplevel(level=0, axis=1)
    summary_savepath = savedir + f'HVTN137_total_protein_pivot_summary.xlsx'
    pivot_summary.to_excel(summary_savepath)


    ## CHECKS ----------------------------------------------------------------##

    # here's BAMA from the same study -- check if same samples were tested
    other = pd.read_csv('/trials/vaccine/p137/s001/qdata/LabData/BAMA_Mucosal/HVTN_137_MUCOSAL_LUM05_20250528.txt', sep="\t")
    assert len(set(outputs.guspec).intersection(other.GUSPEC)) == outputs.guspec.nunique()
    assert len(set(outputs.guspec).symmetric_difference(other.GUSPEC)) == 0

    another = pd.read_csv('/trials/vaccine/p137/s001/qdata/LabData/BAMA_Mucosal/HVTN_137_TOTAL_ANTIBODY_LUM05_20250528.txt', sep="\t")
    assert len(set(outputs.guspec).intersection(another.GUSPEC)) == outputs.guspec.nunique()
    assert len(set(outputs.guspec).symmetric_difference(another.GUSPEC)) == 0

if __name__ == '__main__':
    main()
