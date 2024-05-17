## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 03/21/2024
# Purpose:  - Process 303 NAb Data from doria-rose VRC for stats
#           - Generate pivot table summary
## ---------------------------------------------------------------------------##
import pandas as pd
import datetime, os
import sdmc_tools.process as sdmc
import sdmc_tools.constants as constants

## Read in data --------------------------------------------------------------##
def main():
    yaml_path = os.path.dirname(__file__) + "/paths.yaml"
    with open(yaml_path, 'r') as file:
        yaml_dict = yaml.safe_load(file)

    input_data_path = yaml_dict['input_data_path']
    # input_data_path = "/trials/vaccine/p303/s001/qdata/LabData/VRC_NAb_pass-through/HVTN303_VRC_Doria-Rose_NAb_20240321.xlsx"
    data = pd.read_excel(input_data_path)

    # Reformat input data --------------------------------------------------------##
    # lowercase column names
    data.columns = [i.lower().replace(" ", "_") for i in data.columns]
    renaming = {
        "pseudovirus": "isolate",
        "barcode":"guspec",
        "unit":"result_units",
        "result_id50": "ID50",
        "result_id80":"ID80"
    }
    data = data.rename(columns=renaming)

    # cast titer threshold to long with one result column
    data = data.melt(id_vars=["guspec", "isolate", "assay_date", "result_units", "llod"],
              value_vars=["ID50", "ID80"],
              var_name="titer_threshold",
              value_name="result"
             )

    # Data handler  --------------------------------------------------------------##
    dh = sdmc.DataHandler(
        input_data = data,
        guspec_col = 'guspec',
        network = 'hvtn'
    )

    # merge ldms
    dh.add_ldms(cols = constants.STANDARD_COLS,
                incl_spec_type=True,
                map_drawdt=True,
                relabel=True,
                enforce_typing=True)

    # add human-entered metadata
    dh.add_metadata({
        "network": "HVTN",
        "upload_lab_id": "C7",
        "assay_lab_name": "Doria-Rose Lab",
        "instrument": "Spectramax",
        "assay_type": "Neutralizing Antibody (NAb)",
        "assay_subtype": "TZM-bl",
        "assay_details": "No blocking",
        "specrole": "Sample",
    })

    # add sdmc processing info
    dh.add_sdmc_processing_info(input_data_path=input_data_path)

    # Format final outputs and save ----------------------------------------------##
    output_data = dh.processed_data.copy()
    cols = ['network', 'protocol', 'specrole', 'guspec', 'ptid', 'visitno',
             'drawdt', 'spectype', 'spec_primary', 'spec_additive',
             'spec_derivative', 'upload_lab_id', 'assay_lab_name', 'assay_type',
             'assay_subtype', 'assay_details','instrument', 'isolate', 'assay_date',
             'result_units', 'titer_threshold', 'result', 'llod',
             'sdmc_processing_datetime', 'sdmc_data_receipt_datetime',
             'input_file_name']
    output_data = output_data[cols]

    date = datetime.date.today().isoformat()
    savepath = yaml_dict["savedir"] + f"DRAFT_HVTN303_NAb_VRC_processed_{date}.txt"
    output_data.to_csv(savepath, sep="\t", index=False)

# Generate pivot summary and save --------------------------------------------##
def pivot():
    pivot_summary = pd.pivot_table(
        output_data[["guspec","isolate","titer_threshold", "result"]],
        index="guspec",
        columns=["isolate","titer_threshold"],
        aggfunc='count'
    )
    pivot_summary.columns = pivot_summary.columns.droplevel(level=0)

    pivot_savepath = f"/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN303/assays/nAb_VRC/misc_files/data_processing/HVTN303_NAb_VRC_pivot_summary_{date}.xlsx"
    pivot_summary.to_excel(pivot_savepath)

if __name__ == '__main__':
    main()
    pivot()
