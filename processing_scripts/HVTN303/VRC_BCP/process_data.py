## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 04/19/2024
# Purpose:  - Process HVTN303 B Cell Phenotyping data from Andrews Lab for stats
#           - generate pivot summaries of data
## ---------------------------------------------------------------------------##
import pandas as pd
import datetime, os
import yaml
import sdmc_adhoc_processing.process as sdmc
import sdmc_adhoc_processing.constants as constants
import sdmc_adhoc_processing.utilities as utilities

## custom processing ---------------------------------------------------------##
def main():
    yaml_path = os.path.dirname(__file__) + "/paths.yaml"
    with open(yaml_path, 'r') as file:
        yaml_dict = yaml.safe_load(file)

    input_data_path = yaml_dict['input_data_path']
    # input_data_path = '/trials/vaccine/p303/s001/qdata/LabData/VRC_BCP_pass-through/HVTN303_VRC_Andrews_BCP_20240418.csv'
    input_data = pd.read_csv(input_data_path)

    input_data = input_data.dropna(how='all').dropna(how='all', axis=1)

    # convert wide to long
    data = input_data.melt(id_vars = ['GLOBAL_ID', 'Sample identifier', 'Visit',
                                        'Test Date', 'Test Type'],
                           value_vars = ['B cells counts', 'B cell %', 'IgG count', 'IgG %', 'Trimer4571+ Count',
                                                'Trimer4571+ %', 'Trimer4571+Trimer 6931+ Count',
                                                'Trimer4571+Trimer 6931+ %', 'Trimer4571+Trimer 6931+FP+ Count',
                                                'Trimer4571+Trimer 6931+FP+ %', 'Trimer4571+ FP+ Count',
                                                'Trimer4571+ FP+ Count %', 'FP+ Count', 'FP+ %',
                                                'Trimer 6931+FP+ Count', 'Trimer 6931+FP+ %', 'Trimer 6931+ Count',
                                                'Trimer 6931+ %'],
                           var_name = "subset")

    # standardize names
    data['metric'] = data.subset.str.rpartition(' ', expand=True)[2].map({
        'counts':'result_count',
        'Count':'result_count',
        '%': 'result_percent_from_lab',
        'count':'result_count',
    })
    data.subset = data.subset.str.rpartition(' ', expand=True)[0]

    remap_subset = {
        'B cells': 'B cell',
        'B cell': 'B cell',
        'IgG': 'IgG',
        'Trimer4571+': 'Trimer4571+',
        'Trimer4571+Trimer 6931+': 'Trimer4571+Trimer 6931+',
        'Trimer4571+Trimer 6931+FP+': 'Trimer4571+Trimer 6931+FP+',
        'Trimer4571+ FP+': 'Trimer4571+ FP+',
        'Trimer4571+ FP+ Count': 'Trimer4571+ FP+', #this is a lab typo
        'FP+': 'FP+',
        'Trimer 6931+FP+': 'Trimer 6931+FP+',
        'Trimer 6931+': 'Trimer 6931+'
    }
    data.subset = data.subset.map(remap_subset)

    # standardize column names
    data.columns = [i.lower().replace(" ","_") for i in data.columns]

    # convert to wide on count vs percent
    data = data.pivot(index=['global_id', 'sample_identifier', 'visit', 'test_date', 'test_type', 'subset'],
                     columns=["metric"],
                     values=["value"]
                    ).reset_index()

    # fix column names
    data = data.droplevel(level=1, axis=1)
    data.columns = ['global_id', 'sample_identifier', 'visit', 'test_date', 'test_type', 'subset',
           'result_count', 'result_percent_from_lab']

    # add percent column standardized to percents
    def scientific_to_pct(x):
        x = round(float(x) * 100, 5)
        return str(x) + "%"

    data['result_percent'] = data['result_percent_from_lab']
    data.loc[~data.result_percent.str.contains("%"), "result_percent"] = data.loc[~data.result_percent.str.contains("%")].result_percent.apply(scientific_to_pct)

    data = data.rename(columns={'global_id':'guspec'})
    data = data.drop(columns = ["visit", "sample_identifier"])

    ## standard processing -------------------------------------------------------##
    # hand-entered metadata
    metadata = {"network": "HVTN",
                "upload_lab_id": "C7",
                "assay_lab_name": "Andrews Lab (VRC)",
                "instrument": "BD S6 Sorter",
                "assay_type": "B Cell Phenotyping",
                "specrole": "Sample"}

    #read in ldms
    today = datetime.date.today().strftime("%Y%m%d")
    ldms_path = f"/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN303/specimens/ldms_feed/hvtn.ldms303.{today}.csv"
    ldms = pd.read_csv(ldms_path, usecols=constants.STANDARD_COLS)

    # standard processing
    outputs = sdmc.standard_processing(
        input_data = data,
        input_data_path = input_data_path,
        guspec_col = 'guspec',
        network = "HVTN",
        metadata_dict = metadata,
        ldms = ldms,
        ldms_usecols = constants.STANDARD_COLS
    )

    ## reorder columns and save --------------------------------------------------##
    # manifest = pd.read_csv(
    #     "/networks/vtn/lab/SDMC_labscience/operations/documents/templates/assay/template_testing/512--2023000139.txt",
    #     sep="\t"
    # )
    reorder = ['network', 'protocol', 'specrole', 'guspec', 'ptid', 'visitno',
           'drawdt', 'spectype', 'spec_primary', 'spec_additive',
           'spec_derivative', 'upload_lab_id', 'assay_lab_name', 'assay_type',
           'instrument', 'subset', 'result_count', 'result_percent', 'result_percent_from_lab', 'test_date',
           'test_type', 'sdmc_processing_datetime', 'sdmc_data_receipt_datetime',
           'input_file_name']
    outputs = outputs[reorder]

    date = datetime.date.today().isoformat()
    savepath = savedir + f'DRAFT_HVTN303_VRC_BCP_processed_{date}.txt'
    outputs.to_csv(savepath, index=False, sep="\t")

## output pivot summary ------------------------------------------------------##
def pivot():
    summary = outputs.melt(id_vars=['guspec','subset'], value_vars=['result_count','result_percent'], var_name="result_units")
    summary = summary.pivot_table(values='value', index='guspec', columns=['subset','result_units'], aggfunc='count')
    summary.to_excel("/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN303/assays/BCP/misc_files/data_processing/HVTN303_VRC_BCP_summary.xlsx")

if __name__ == '__main__':
    main()
    pivot()
