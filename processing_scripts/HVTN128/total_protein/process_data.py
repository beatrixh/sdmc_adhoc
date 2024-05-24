## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 04/09/2024
# Purpose:  - Process HVTN128 total protein data from tomaras lab
#           - Generate pivot summary
## ---------------------------------------------------------------------------##
import pandas as pd
import datetime, os
import yaml
import sdmc_tools.process as sdmc
import sdmc_tools.constants as constants
## ---------------------------------------------------------------------------##
def main():
    yaml_path = "/home/bhaddock/repos/sdmc-adhoc/processing_scripts/HVTN128/total_protein" + "/paths.yaml"
    with open(yaml_path, 'r') as file:
        yaml_dict = yaml.safe_load(file)

    input_data_path = yaml_dict['input_data_path']
    total_protein = pd.read_csv(input_data_path)

    total_protein = total_protein.rename(columns={
        "global_identifier":"guspec",
        "observed_concentration":"concentration"
        })

    ldms = pd.read_csv(constants.LDMS_PATH_HVTN, usecols=constants.STANDARD_COLS)
    ldms = ldms.loc[ldms.lstudy==128.].drop_duplicates()
    ldms = ldms.loc[ldms.guspec.isin(total_protein.guspec.unique().tolist())]

    metadata_dict = {
        "network": "HVTN",
        "upload_lab_id": "GT",
        "assay_lab_name": "Tomaras Lab (Duke)",
        "assay_kit": "Quant-iT",
        "assay_type": "Total Protein",
        "specrole": "Sample",
    }

    outputs = sdmc.standard_processing(
        input_data=total_protein,
        input_data_path=input_data_path,
        guspec_col="guspec",
        network='hvtn',
        metadata_dict=metadata_dict,
        ldms=ldms
    )

    if len(outputs.spectype.isna().sum()) > 0:
        raise Exception("THERE ARE NANS IN SPEDCTYPE; NEED TO FIX")

    outputs = outputs.drop(columns=[
        "assay_platform", "lab", "study", "visit_identifier",
        "timepoint", "timepoint_units", "isotype", "sample_identifier"
    ])

    today = datetime.date.today().isoformat()
    savedir = yaml_dict['savedir']
    savepath = savedir + f'DRAFT_HVTN128_total_protein_processed_{today}.txt'
    outputs.to_csv(savepath, index=False, sep="\t")
## ---------------------------------------------------------------------------##
def pivot():
    pivot_summary = pd.pivot_table(outputs[["ptid","visitno","spectype","concentration"]],
                                   index=["ptid", "visitno"],
                                   columns="spectype",
                                   aggfunc='count',
                                   dropna=False,
                                   fill_value=0
                                  )

    pivot_summary = pivot_summary.droplevel(level=0, axis=1)
    summary_savepath = savedir + f'HVTN128_total_protein_pivot_summary.xlsx'
    pivot_summary.to_excel(summary_savepath)

    ## Checks --------------------------------------------------------------------##
    manifest = pd.read_excel("/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN128/assays/smc/manifests/Precision_HVTN_128_GT_28SEP22_1342.xlsx")
    manifest["guspec"] = manifest["Original Id"]
    manifest.columns = [i.lower().replace(" ","_") for i in manifest.columns]

    # counts of material types of things in the manifest but NOT in the smc data
    print("counts of material types of rows in the manifest but NOT in the smc data")
    print(manifest.loc[~(manifest.guspec.isin(outputs.guspec))].material_type.value_counts())

    print("counts of material types of things from the manifest IN the smc data")
    print(manifest.loc[(manifest.guspec.isin(outputs.guspec))].material_type.value_counts())

    print("guspec of the rectal biopsy guspec thats in the manifest but NOT in the smc data:")
    print(manifest.loc[~(manifest.guspec.isin(outputs.guspec)) & (manifest.material_type=="Rectal Biopsy")].guspec)

    if len(set(manifest.subject_id.astype(int)).symmetric_difference(outputs.ptid.astype(int))) == 0:
        print("we have a match between prtids in manifest and results")
    else:
        print("mismatch between results and manifest")

if __name__ == '__main__':
    main()
    pivot()
