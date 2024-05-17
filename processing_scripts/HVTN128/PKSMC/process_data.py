## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 04/09/2024
# Purpose:  - Process HVTN128 pk-smc data from tomaras lab
#           - Generate pivot summary
## ---------------------------------------------------------------------------##
import pandas as pd
import datetime, os
import yaml
import sdmc_tools.process as sdmc
import sdmc_tools.constants as constants
## ---------------------------------------------------------------------------##
def main():
    yaml_path = os.path.dirname(__file__) + "/paths.yaml"
    with open(yaml_path, 'r') as file:
        yaml_dir = yaml.safe_load(file)
    # new_pksmc_input_path = '/trials/vaccine/p128/s001/qdata/LabData/PK-SMC_pass-through/HVTN128_PK-SMC_20240402_A.csv'
    new_pksmc_input_path = yaml_dir["input_data_path"]
    new_pksmc = pd.read_csv(new_pksmc_input_path)

    new_pksmc.columns = [i.lower().replace(" ", "_") for i in new_pksmc.columns]

    new_pksmc = new_pksmc.rename(columns = {
        "global_identifier": "guspec",
        "observed_concentration":"concentration",
        "lower_limit_of_detection": "llod",
        "lower_limit_of_quantitation":"lloq",
        "lower_limit_of_detection_units": "llod_units",
        "lower_limit_of_quantitation_units": "lloq_units"
    })

    new_pksmc = new_pksmc.drop(columns = ["lab", "study",
                                          "visit_identifier",
                                          "sample_identifier", "timepoint",
                                          "timepoint_units", "timepoint_baseline",
                                          "assay_platform", "analysis_software"
                                         ])
    ## ---------------------------------------------------------------------------##
    ldms = pd.read_csv(constants.LDMS_PATH_HVTN, usecols=constants.STANDARD_COLS)
    ldms = ldms.loc[ldms.lstudy==128.].drop_duplicates()
    dh = sdmc.DataHandler(
        input_data = new_pksmc,
        guspec_col = 'guspec',
        network = 'hvtn'
    )

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

    dh.ldms = ldms.copy()

    # merge ldms
    dh.add_ldms(cols = constants.STANDARD_COLS,
                incl_spec_type=True,
                map_drawdt=True,
                relabel=True,
                enforce_typing=True)
    dh.processed_data.spectype = dh.processed_data.apply(
        lambda x: constants.SPEC_TYPE_DEFN_MAP[x.spec_primary, x.spec_derivative],
        axis=1
    )
    if dh.processed_data.spectype.isna().sum() > 0:
        raise Exception("THERE ARE NANS IN SPEDCTYPE; NEED TO FIX")

    dh.add_metadata({
        "network": "HVTN",
        "upload_lab_id": "GT",
        "assay_lab_name": "Tomaras Lab (Duke)",
        "instrument": "SMCxPRO",
        "lab_software_version": "xPRO 1.2.69",
        "assay_type": "PK-SMC",
        "specrole": "Sample",
    })

    dh.add_sdmc_processing_info(input_data_path=new_pksmc_input_path)
    ## ---------------------------------------------------------------------------##
    outputs = dh.processed_data.copy()
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
        'instrument',
        'lab_software_version',
        'analyte',
        'analyte_lot',
        'concentration_method',
        'concentration_standard',
        'concentration',
        'concentration_units',
        'curve_location',
        'dataset_creation_date',
        'dataset_version',
        'description',
        'dilution',
        'exclude',
        'expected_concentration',
        'experiment_date',
        'experiment_identifier',
        'llod',
        'llod_units',
        'lloq',
        'lloq_units',
        'plate_identifier',
        'plate_type',
        'response',
        'response_coefficient_of_variation',
        'response_standard_deviation',
        'run_identifier',
        'sample_contamination',
        'sample_role',
        'sample_species',
        'sample_type',
        'source_file_short',
        'summary',
        'summary_method',
        'summary_type',
        'well',
        'sdmc_processing_datetime',
        'sdmc_data_receipt_datetime',
        'input_file_name'
    ]

    if set(reorder).symmetric_difference(outputs.columns):
        raise Exception("TRYING TO DROP OR ADD COLUMNS ON REORDER; CHECK")

    outputs = outputs[reorder]

    ## save data -----------------------------------------------------------------##
    today = datetime.date.today().isoformat()
    savepath = yaml_dir["savedir"] + f'DRAFT_HVTN128_PKSMC_processed_{today}.txt'
    outputs.to_csv(savepath, index=False, sep="\t")

## pivot summary -------------------------------------------------------------##
def pivot():
    pivot_summary = pd.pivot_table(outputs[["ptid","visitno","spectype","concentration"]],
                                   index=["ptid", "visitno"],
                                   columns="spectype",
                                   aggfunc='count',
                                   dropna=False,
                                   fill_value=0
                                  )
    pivot_summary = pivot_summary.droplevel(level=0, axis=1)

    summary_savepath = f'/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN128/assays/smc/misc_files/data_processing/HVTN128_pksmc_pivot_summary.xlsx'
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
