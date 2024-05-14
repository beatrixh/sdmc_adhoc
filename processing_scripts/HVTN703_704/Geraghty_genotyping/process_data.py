## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 02/29/2024
# Purpose:  - Process HLA genotyping data for 704 from Geraghty lab for stats
#           - generate pivot summaries of data
## ---------------------------------------------------------------------------##
import pandas as pd
import datetime, os
import yaml
import sdmc_tools.process as sdmc
import sdmc_tools.constants as constants

## Read in data --------------------------------------------------------------##
def main():
    yaml_path = os.path.dirname(__file__) + "/paths.yaml"
    with open(yaml_path, 'r') as file:
        yaml_dict = yaml.safe_load(file)

    input_data_path = yaml_dict['input_data_path']
    # input_data_path = '/trials/vaccine/p704/s670/qdata/HLA_pass-through/HVTN704-HLA-CCR5.xlsx'
    hla_sheet = pd.read_excel(input_data_path, sheet_name="HLA")
    ccr5_sheet = pd.read_excel(input_data_path, sheet_name="CCR5")

    bmetadata = ccr5_sheet.iloc[22:]
    ccr5_sheet = ccr5_sheet.iloc[:18]

    ## Reshape/format input data -------------------------------------------------##

    # add guspec
    hla_sheet["guspec"] = hla_sheet["Sample ID"].apply(lambda x: "-".join(x.split("-")[-2:]))
    ccr5_sheet["guspec"] = ccr5_sheet["Sample ID"].apply(lambda x: "-".join(x.split("-")[-2:]))

    # convert to lower_case
    hla_sheet.columns = [i.lower().replace(" ", "_") for i in hla_sheet.columns]
    ccr5_sheet.columns = [i.lower().replace(" ", "_") for i in ccr5_sheet.columns]

    # rename cols
    hla_sheet_col_map = {
        'diploid_ambiguities': 'ambiguities_diploid',
        'allele_1_ambiguities': 'ambiguities_allele_1',
        'allele_2_ambiguities': 'ambiguities_allele_2'
    }

    hla_sheet = hla_sheet.rename(columns = hla_sheet_col_map)

    ccr5_sheet_col_map = {
        "allele1": "allele_1",
        "allele2": "allele_2"
    }

    ccr5_sheet = ccr5_sheet.rename(columns = ccr5_sheet_col_map)

    # add clarifying columns to ccr5_sheet
    ccr5_sheet_remapping = {
        "normal": "no deletion of 32bp",
        "deletion": "deletion of 32bp (CCR5-delta 32)"
    }

    ccr5_sheet["allele_1_detail"] = ccr5_sheet.allele_1.map(ccr5_sheet_remapping)
    ccr5_sheet["allele_2_detail"] = ccr5_sheet.allele_2.map(ccr5_sheet_remapping)

    # add subtypes
    hla_sheet["assay_subtype"] = "HLA"
    ccr5_sheet["assay_subtype"] = "CCR5"
    ## Merge on human-entered metadata -------------------------------------------##

    metadata = pd.DataFrame({
        "network": ["HVTN"],
        "upload_lab_id": ["DG"],
        "assay_lab_name": ["Geraghty Lab (FHCRC)"],
        'instrument': ["Illumina NGS"],
        'assay_type': ["Genotyping"],
        'specrole': ["Sample"],
    })

    hla_sheet = hla_sheet.merge(metadata, how = 'cross')
    ccr5_sheet = ccr5_sheet.merge(metadata, how = 'cross')

    ## Merge on ldms -------------------------------------------------------------##

    ldms_usecols = ["lstudy",
                    "addstr",
                    "dervstr",
                    "primstr",
                    "txtpid",
                    "vidval",
                    "drawdm", "drawdd", "drawdy", "guspec"]
    ldms_col_map = {
        "lstudy": "protocol",
        "addstr": "spec_additive",
        "dervstr": "spec_derivative",
        "primstr": "spec_primary",
        "txtpid": "ptid",
        "vidval": "visitno",
    }
    spec_map = {
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
        ('VCS', 'SWB'): 'Cervicovaginal Secretions (Swab)'
    }

    # read in ldms
    hvtn_csdb_path = '/data/pipelines/fstrf-ldms-postimport/hvtn/hvtn.csdb.csv'
    ldms = pd.read_csv(hvtn_csdb_path, usecols=ldms_usecols)

    #subset down
    gus = list(hla_sheet.guspec) + list(ccr5_sheet.guspec)
    ldms = ldms.loc[ldms.guspec.isin(gus)]

    # reformat
    ldms = ldms.rename(columns=ldms_col_map)

    # add date column
    ldms["drawdt"] = ldms.apply(
        lambda x: datetime.date(x.drawdy, x.drawdm, x.drawdd).isoformat(), axis=1
    )
    ldms = ldms.drop(columns=["drawdy", "drawdm", "drawdd"])

    def map_spectype(x):
        try:
            return spec_map[x.spec_primary, x.spec_derivative]
        except:
            print(f"{x} missing from spec map!")
            return "MISSING FROM MAP"

    ldms.loc[ldms.spec_derivative.isna(), "spec_derivative"] = "N/A"
    ldms["spectype"] = ldms.apply(lambda x: map_spectype(x), axis=1)

    ldms.ptid = ldms.ptid.astype(int).astype(str)
    ldms.protocol = ldms.protocol.astype(int).astype(str)

    hla_sheet = hla_sheet.merge(ldms, on="guspec", how="left")
    ccr5_sheet = ccr5_sheet.merge(ldms, on="guspec", how="left")

    ## Add SDMC processing info --------------------------------------------------##

    sdmc_processing_datetime = datetime.datetime.now().replace(microsecond=0).isoformat()
    data_receipt_datetime = datetime.datetime.fromtimestamp(os.path.getmtime(input_data_path)).replace(microsecond=0).isoformat()

    sdmc_metadata = pd.DataFrame({
        "sdmc_processing_datetime": [sdmc_processing_datetime],
        # "sdmc_processing_version": [1.0],
        "sdmc_data_receipt_datetime": [data_receipt_datetime],
        "input_file_name": [input_data_path.split("/")[-1]],
    })

    hla_sheet = hla_sheet.merge(sdmc_metadata, how="cross")
    ccr5_sheet = ccr5_sheet.merge(sdmc_metadata, how="cross")

    ## Reformat and save ---------------------------------------------------------##

    acols = ['network', 'protocol', 'specrole', 'guspec', 'ptid', 'visitno',
             'drawdt', 'spectype', 'spec_primary', 'spec_additive',
             'spec_derivative', 'upload_lab_id', 'assay_lab_name', 'assay_type',
             'assay_subtype', 'instrument',
             'sample_id', 'locus', 'allele_1', 'allele_2', 'comments',
             'ambiguities_diploid', 'ambiguities_allele_1', 'ambiguities_allele_2',
             'sdmc_processing_datetime', 'sdmc_data_receipt_datetime',
             'input_file_name']

    bcols = ['network', 'protocol', 'specrole', 'guspec', 'ptid', 'visitno',
             'drawdt', 'spectype', 'spec_primary', 'spec_additive',
             'spec_derivative', 'upload_lab_id', 'assay_lab_name', 'assay_type',
             'assay_subtype', 'instrument',
             'sample_id', 'allele_1', 'allele_2', 'allele_1_detail', 'allele_2_detail',
             'sdmc_processing_datetime', 'sdmc_data_receipt_datetime',
             'input_file_name']

    acol_mismatch = set(acols).symmetric_difference(hla_sheet.columns)
    bcol_mismatch = set(bcols).symmetric_difference(ccr5_sheet.columns)
    if len(acol_mismatch) > 0:
        raise Exception(f"HLA colummn mismatch: {acol_mismatch}")
    if len(bcol_mismatch) > 0:
        raise Exception(f"CCR5 colummn mismatch: {bcol_mismatch}")

    hla_sheet = hla_sheet[acols]
    ccr5_sheet = ccr5_sheet[bcols]

    timestamp = datetime.datetime.today().strftime('%Y-%m-%d')
    # output_dir = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN703_704/assays/genotyping/misc_files/data_processing/'
    output_dir = yaml_dict['savedir']

    fname1 = f"DRAFT_HVTN704-HLA_processed_{timestamp}.txt"
    hla_sheet.to_csv(output_dir + fname1, sep="\t", index=False)

    fname2 = f"DRAFT_HVTN704-CCR5_processed_{timestamp}.txt"
    ccr5_sheet.to_csv(output_dir + fname2, sep="\t", index=False)

## Output pivot summaries ----------------------------------------------------##

def pivot():
    hla_sheet_summary = hla_sheet.copy()
    hla_sheet_summary = hla_sheet_summary.melt(id_vars = ['network', 'protocol', 'specrole', 'guspec', 'ptid', 'visitno',
           'drawdt', 'spectype', 'spec_primary', 'spec_additive',
           'spec_derivative', 'upload_lab_id', 'assay_lab_name', 'assay_type',
           'instrument', 'sample_id', 'locus', 'comments',
           'ambiguities_diploid', 'ambiguities_allele_1', 'ambiguities_allele_2',
           'sdmc_processing_datetime', 'sdmc_data_receipt_datetime',
           'input_file_name'], value_vars=["allele_1", "allele_2"], var_name="allele", value_name = "value")
    hla_sheet_summary = hla_sheet_summary.pivot_table(index="guspec", columns=["locus", "allele"], aggfunc = ("count"))[["value"]]
    hla_sheet_summary = hla_sheet_summary.droplevel(level=0, axis=1)
    hla_sheet_summary.to_excel(f"{output_dir}HLA_dataset_summary_{timestamp}.xlsx", sheet_name="HLA")

    ccr5_sheet_summary = ccr5_sheet.groupby(["guspec"]).count()[["allele_1", "allele_2"]]
    ccr5_sheet_summary.to_excel(f"{output_dir}CCR5_dataset_summary_{timestamp}.xlsx", sheet_name="CCR5")

if __name__ == '__main__':
    main()
    pivot()
