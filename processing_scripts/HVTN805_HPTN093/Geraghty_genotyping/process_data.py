## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 04/08/2024
# Purpose:  - Process new and old HLA genotyping data for 805 from
#             Geraghty lab for stats
## ---------------------------------------------------------------------------##
import pandas as pd
from datetime import date
import sdmc_adhoc_processing.process as sdmc
import sdmc_adhoc_processing.constants as constants

## ---------------------------------------------------------------------------##
# We need to process both the new and the old data in the same manner, but append
# different "data reciept" dates. The data is formatted the same, so wrapping a
# function below to apply to both new and old data.
def main():
    def process_hla_data(input_data_path, ldms):
        hla = pd.read_excel(input_data_path, sheet_name = "HLA")
        hla.columns = [i.lower().replace(" ","_") for i in hla.columns]

        hla = hla.fillna("")

        hla_rename = {
            'diploid_ambiguities':'ambiguities_diploid',
            'allele_1_ambiguities':'ambiguities_allele_1',
            'allele_2_ambiguities':'ambiguities_allele_2',
        }

        hla = hla.rename(columns = hla_rename)

        # load data handler
        dh_hla = sdmc.DataHandler(
            input_data = hla,
            guspec_col = 'sample_id',
            network = 'hvtn'
        )

        # load ldms
        dh_hla.ldms = ldms

        # merge ldms
        dh_hla.add_ldms(cols = constants.STANDARD_COLS,
                    incl_spec_type=True,
                    map_drawdt=True,
                    relabel=True)

        # add human-entered metadata
        dh_hla.add_metadata({
            "network": "HVTN",
            "upload_lab_id": "DG",
            "assay_lab_name": "Geraghty Lab (FHCRC)",
            "instrument": "Illumina NGS",
            "assay_type": "Genotyping",
            "assay_subtype": "HLA",
            "specrole": "Sample",
        })

        dh_hla.add_sdmc_processing_info(input_data_path=input_data_path)
        hlacols = ['network', 'protocol', 'specrole', 'guspec', 'ptid', 'visitno',
                 'drawdt', 'spectype', 'spec_primary', 'spec_additive',
                 'spec_derivative', 'upload_lab_id', 'assay_lab_name', 'assay_type', 'assay_subtype', 'instrument',
                 'sample_id', 'locus', 'allele_1', 'allele_2', 'comments',
                 'ambiguities_diploid', 'ambiguities_allele_1', 'ambiguities_allele_2',
                 'sdmc_processing_datetime', 'sdmc_data_receipt_datetime',
                 'input_file_name']

        return dh_hla.processed_data[hlacols].drop_duplicates()

    def process_ccr5_data(input_data_path, ccr5, ldms):
        ccr5.columns = [i.lower().replace(" ","_") for i in ccr5.columns]
        ccr5_rename = {
            'allele1': 'allele_1',
            'allele2': 'allele_2'
        }

        ccr5 = ccr5.rename(columns = ccr5_rename)
        notation_map = {
            "normal": "no deletion of 32bp",
            "deletion": "deletion of 32bp (CCR5-delta 32)"
        }
        ccr5["allele_1_detail"] = ccr5.allele_1.map(notation_map)
        ccr5["allele_2_detail"] = ccr5.allele_2.map(notation_map)

        dh_ccr = sdmc.DataHandler(
            input_data = ccr5,
            guspec_col = 'sample_id',
            network = 'hvtn'
        )

        # load ldms
        # dh.load_ldms(usecols=constants.STANDARD_COLS)
        dh_ccr.ldms = ldms

        # merge ldms
        dh_ccr.add_ldms(cols = constants.STANDARD_COLS,
                    incl_spec_type=True,
                    map_drawdt=True,
                    relabel=True)

        # add human-entered metadata
        dh_ccr.add_metadata({
            "network": "HVTN",
            "upload_lab_id": "DG",
            "assay_lab_name": "Geraghty Lab (FHCRC)",
            "instrument": "Illumina NGS",
            "assay_type": "Genotyping",
            "assay_subtype": "CCR5",
            "specrole": "Sample",
        })

        dh_ccr.add_sdmc_processing_info(input_data_path=input_data_path)

        ccr5_cols = ['network', 'protocol', 'specrole', 'guspec', 'ptid', 'visitno',
                 'drawdt', 'spectype', 'spec_primary', 'spec_additive',
                 'spec_derivative', 'upload_lab_id', 'assay_lab_name', 'assay_type', 'assay_subtype', 'instrument',
                 'allele_1', 'allele_2', 'allele_1_detail','allele_2_detail',
                 'sdmc_processing_datetime', 'sdmc_data_receipt_datetime',
                 'input_file_name']

        return dh_ccr.processed_data[ccr5_cols].drop_duplicates()
    ## read in data --------------------------------------------------------------##
    input_data_path_new = '/trials/vaccine/p805/s001/qdata/LabData/Genotyping_pass-through/HVTN805-20240405-additional-gID.xlsx'
    input_data_path_old = "/trials/vaccine/p805/s001/qdata/LabData/Genotyping_pass-through/HVTN805-HLA-CCR5-name-revised.xlsx"

    ldms = pd.read_csv(constants.LDMS_PATH_HVTN)
    ldms = ldms.loc[ldms.lstudy==805.]
    ldms = ldms.drop_duplicates()

    ## process data --------------------------------------------------------------##
    old_hla = process_hla_data(input_data_path_old, ldms)
    new_hla = process_hla_data(input_data_path_new, ldms)
    hla = pd.concat([old_hla, new_hla])

    # ccr5 data has extraneous values that have to be trimmed
    # (they're merged on as "detail" columns)
    ccr5_old = pd.read_excel(input_data_path_old, sheet_name="CCR5")
    ccr5_old = ccr5_old.iloc[:10]

    ccr5_new = pd.read_excel(input_data_path_new, sheet_name="CCR5")
    ccr5_new = ccr5_new.iloc[:1,:3]

    ccr5_old = process_ccr5_data(input_data_path_old, ccr5_old, ldms)
    ccr5_new = process_ccr5_data(input_data_path_new, ccr5_new, ldms)
    ccr5 = pd.concat([ccr5_old, ccr5_new])

    ## save to tab-delimited txt -------------------------------------------------##
    savedir = "/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN805_HPTN093/assays/genotyping/misc_files/data_processing/"
    today = datetime.date.today().isoformat()

    hla.to_csv(savedir + f"DRAFT_HVTN805_HLA_genotyping_processed_{today}.txt",
               sep="\t",
               index=False
              )
    ccr5.to_csv(savedir + f"DRAFT_HVTN805_CCR5_genotyping_processed_{today}.txt",
               sep="\t",
               index=False
              )
    ## checks --------------------------------------------------------------------##
    if set(hla.guspec).symmetric_difference(ccr5.guspec):
        raise Exception("guspecs dont match between hla and ccr5")

    print("guspec value counts for hla:")
    print(hla.guspec.value_counts())

    print("guspec value counts for ccr5:")
    print(ccr5.guspec.value_counts())

## pivot summary -------------------------------------------------------------##
def pivot():
    hla_summary = hla.melt(id_vars=["guspec","locus"],
                           value_vars=["allele_1", "allele_2"],
                           var_name="allele")
    hla_summary= pd.pivot_table(hla_summary,
                                values="value",
                                index="guspec",
                                columns=["locus", "allele"],
                                aggfunc='count')
    hla_summary.to_excel(savedir + f"HLA_pivot_summary_{today}.xlsx")

    ccr5_summary = ccr5.melt(id_vars=["guspec"],
                             value_vars=["allele_1", "allele_2"],
                             var_name="allele")
    ccr5_summary = pd.pivot_table(ccr5_summary,
                                  values="value",
                                  index="guspec",
                                  columns="allele",
                                  aggfunc='count')
    ccr5_summary.to_excel(savedir + f"CCR5_pivot_summary_{today}.xlsx")

if __name__ == '__main__':
    main()
    pivot()
