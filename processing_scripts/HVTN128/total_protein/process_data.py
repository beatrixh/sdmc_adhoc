## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 04/09/2024
# Purpose:  - Process HVTN128 total protein data from tomaras lab
#           - Generate pivot summary
## ---------------------------------------------------------------------------##
import pandas as pd
import datetime, os
import yaml
import sdmc_adhoc_processing.process as sdmc
import sdmc_adhoc_processing.constants as constants
## ---------------------------------------------------------------------------##
def main():
    yaml_path = os.path.dirname(__file__) + "/paths.yaml"
    with open(yaml_path, 'r') as file:
        yaml_dict = yaml.safe_load(file)

    input_data_path = yaml_dict['input_data_path']
    total_protein.columns = [i.lower().replace(" ","_") for i in total_protein.columns]

    total_protein = total_protein.rename(columns={
        "global_identifier":"guspec",
        "observed_concentration":"concentration"
        })

    ldms = pd.read_csv(constants.LDMS_PATH_HVTN, usecols=constants.STANDARD_COLS)
    ldms = ldms.loc[ldms.lstudy==128.].drop_duplicates()

    dh = sdmc.DataHandler(
        input_data = total_protein,
        guspec_col = 'guspec',
        network = 'hvtn'
    )
    dh.ldms = ldms.copy()

    # merge ldms
    dh.add_ldms(cols = constants.STANDARD_COLS,
                incl_spec_type=True,
                map_drawdt=True,
                relabel=True,
                enforce_typing=True)

    constants.SPEC_TYPE_DEFN_MAP = {
        ('BLD','PLA'): 'Plasma',
        ('BLD','SER'): 'Serum',
        ('REC','SPG'): 'Rectal Sponge',
        ('REC','SUP'): 'Rectal Biopsy',
        ('CER','SPG'): 'Cervical Sponge',
        ('CER','SUP'): 'Cervical Biopsy',
        ('VAG','SUP'): 'Vaginal Biopsy',
        ('VAG','WCK'): 'Vaginal Weck',
        ('SAL','FLD'): 'Saliva',
        ('SAL','SAL'): 'Saliva',
        ('SEM','SEM'): 'Semen',
        ('SEM','FLD'): 'Semen',
        ('BLD','CEL'): 'PBMC',
        ('BLD', 'CSR'): 'Serum',
        ('BLD', 'DBS'): 'Dried Blood Spot',
        ('BLD', 'LYS'): 'Whole Blood (Lysed)',
        ('VCS', 'FLD'): 'Cervicovaginal Secretions (Fluid)',
        ('VCS', 'MUC'): 'Cervicovaginal Secretions (Mucus)',
        ('VAG', 'SWB'): 'Vaginal Swab',
        ('BLD', 'BLD'): 'Whole Blood',
        ('VCS', "N/A"): 'Cervicovaginal Secretions',
        ('BLD', "N/A"): 'Whole Blood',
        ('VCS', 'SWB'): 'Cervicovaginal Secretions (Swab)',
        ('REC', 'BPS'): 'Rectal Biopsy',
        ('VAG', 'BPS'): 'Vaginal Biopsy',
        ('CER', 'BPS'): 'Cervical Biopsy',
        ('REC', 'SEC'): 'Rectal Secretions'
    }
    dh.processed_data.spectype = dh.processed_data.apply(
        lambda x: constants.SPEC_TYPE_DEFN_MAP[x.spec_primary, x.spec_derivative],
        axis=1
    )

    if len(dh.processed_data.spectype.isna().sum()) > 0:
        raise Exception("THERE ARE NANS IN SPEDCTYPE; NEED TO FIX")

    dh.add_metadata({
        "network": "HVTN",
        "upload_lab_id": "GT",
        "assay_lab_name": "Tomaras Lab (Duke)",
        "assay_kit": "Quant-iT",
        "assay_type": "Total Protein",
        "specrole": "Sample",
    })

    dh.add_sdmc_processing_info(input_data_path=input_data_path)

    outputs = dh.processed_data.copy()

    outputs = outputs.drop(columns=[
        "assay_platform", "lab", "study", "visit_identifier",
        "timepoint", "timepoint_units", "isotype", "sample_identifier"
    ])

    reorder = [
        'network',
        'protocol',
        'specrole',
        'guspec',
        'ptid',
        'visitno',
        'drawdt',
        'spectype',
        'spec_primary',
        'spec_additive',
        'spec_derivative',
        'upload_lab_id',
        'assay_lab_name',
        'assay_type',
        'assay_kit',
        'analyte',
        'concentration_units',
        'description',
        'experiment_date',
        'experiment_identifier',
        'concentration',
        'sample_role',
        'sample_type',
        'sdmc_processing_datetime',
        'sdmc_data_receipt_datetime',
        'input_file_name'
    ]
    if set(reorder).symmetric_difference(outputs.columns):
        raise Exception("TRYING TO DROP OR ADD COLUMNS ON REORDER; CHECK")

    outputs = outputs[reorder]

    today = datetime.date.today().isoformat()
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
    summary_savepath = f'/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN128/assays/total_protein/misc_files/data_processing/HVTN128_total_protein_pivot_summary.xlsx'
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

if __name__ == '__main__':
    main()
    pivot()
