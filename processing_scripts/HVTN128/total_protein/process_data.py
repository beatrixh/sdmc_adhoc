## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 09/16/2024
# Purpose:
## ---------------------------------------------------------------------------##
import pandas as pd
import datetime, os
import yaml
import sdmc_tools.refactor_process as sdmc
import sdmc_tools.constants as constants
## ---------------------------------------------------------------------------##

def main():
    input_data_path = "/trials/vaccine/p128/s001/qdata/LabData/Total_Protein_pass-through/HVTN128_Total_Protein_20240402_A.csv"
    original = pd.read_csv(input_data_path)

    new_input_data_path = "/trials/vaccine/p128/s001/qdata/LabData/Total_Protein_pass-through/HVTN128_Total_Protein_SER_20240912_A.csv"
    new = pd.read_csv(new_input_data_path)

    mismatch = set(original.columns).symmetric_difference(new.columns)
    if len(mismatch) > 0:
        raise Exception("Old and new datasets have differing columns")

    ## ---------------------------------------------------------------------------##
    ldms = pd.read_csv(constants.LDMS_PATH_HVTN, usecols=constants.STANDARD_COLS, dtype=constants.LDMS_DTYPE_MAP)
    ldms = ldms.loc[ldms.lstudy==128.].drop_duplicates()

    ## ---------------------------------------------------------------------------##
    def process(data, input_data_path):
        data = data.rename(columns={
        "global_identifier":"guspec",
        "observed_concentration":"concentration"
        })

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
        return outputs

    og_outputs = process(original, input_data_path)
    new_outputs = process(new, new_input_data_path)
    outputs = pd.concat([og_outputs, new_outputs])

    outputs = outputs.drop(columns=[
        "assay_platform", "lab", "study", "visit_identifier",
        "timepoint", "timepoint_units", "isotype", "sample_identifier"
    ])

    ## ---------------------------------------------------------------------------##
    savedir = "/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN128/assays/total_protein/misc_files/data_processing/"
    today = datetime.date.today().isoformat()

    savepath = savedir + f'DRAFT_HVTN128_total_protein_processed_{today}.txt'
    outputs.to_csv(savepath, index=False, sep="\t")

    ## ---------------------------------------------------------------------------##
def pivot():
    summary = pd.pivot_table(outputs,
                             index=['ptid','visitno'],
                             columns='spectype',
                             aggfunc='count',
                             dropna=False,
                             fill_value=0
                            )['concentration']
    summary.to_excel(savedir + "HVTN128_total_protein_sample_summary.xlsx")

    summary2 = pd.pivot_table(outputs,
                             index=['ptid','visitno'],
                             columns=['input_file_name','spectype'],
                             aggfunc='count',
                             dropna=False,
                             fill_value=0
                            )['concentration']
    summary2.to_excel(savedir + "HVTN128_total_protein_sample_summary_by_dataset.xlsx")

    ## Checks --------------------------------------------------------------------##
def checks():
    manifest = pd.read_excel("/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN128/assays/smc/manifests/Precision_HVTN_128_GT_28SEP22_1342.xlsx")
    manifest["guspec"] = manifest["Original Id"]
    manifest.columns = [i.lower().replace(" ","_") for i in manifest.columns]

    ## from new data
    serum_mismatch= set(new_outputs.guspec).symmetric_difference(manifest.loc[manifest.material_type=="Serum"].guspec)
    if len(serum_mismatch) > 0:
        raise Exception("Serum samples in data don't match those in manifest")

    # original checks
    # counts of material types of things in the manifest but NOT in the smc data
    print("counts of material types of rows in the manifest but NOT in the smc data")
    print(manifest.loc[~(manifest.guspec.isin(og_outputs.guspec))].material_type.value_counts())

    print("counts of material types of things from the manifest IN the smc data")
    print(manifest.loc[(manifest.guspec.isin(og_outputs.guspec))].material_type.value_counts())

    print("guspec of the rectal biopsy guspec thats in the manifest but NOT in the smc data:")
    print(manifest.loc[~(manifest.guspec.isin(og_outputs.guspec)) & (manifest.material_type=="Rectal Biopsy")].guspec)

    if len(set(manifest.subject_id.astype(int)).symmetric_difference(og_outputs.ptid.astype(int))) == 0:
        print("we have a match between prtids in manifest and results")
    else:
        print("mismatch between results and manifest")


if __name__=="__main__":
    main()
    pivot()
    checks()
