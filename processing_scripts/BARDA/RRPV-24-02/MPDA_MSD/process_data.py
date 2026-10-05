import pandas as pd
import numpy as np
import os
from datetime import date

import sdmc_tools.process as sdmc_tools
import sdmc_tools.access_ldms as access_ldms

path2 = '/networks/vtn/lab/SDMC_labscience/assays/Binding_Ab_via_MSD/PPD_ThermoFisher/example_data/RRPV-24-02_PPDVAC167_UBT_20260921_v1.csv'
ubt = pd.read_csv(path2)

col_rename = {
    'PROTOCOL': 'protocol',
    'MATRIX': 'spectype',
    'ANALYTE': 'analyte',
    'BARCODE': 'sample_id_lab',
    'SUBJECT': 'ptid',
    'SITE': 'site_id',
    'VISIT': 'visitno',
    'COLLECT_DATE': 'drawdt',
    'LLOQ': 'lloq',
    'Analytical ULOQ': 'uloq_analytical',
    'Dynamic ULOQ': 'uloq_dynamic',
    'RESULT': 'result',
    'RESULT_UNIT': 'result_units',
    'COMMENTS': 'comments_lab',
}

data = pd.concat([recont, ubt]).rename(columns=col_rename)
data = ubt.copy().rename(columns=col_rename)

data.head()

md = {
    'network':'BARDA',
    'specrole':'Sample',
    'upload_lab_id':'P8',
    'assay_lab_name':'PPD-BioA',
    'assay_type':'Binding Antibody',
    'assay_subtype':'MSD',
    'instrument':'APRIL TO ASK ABOUT THIS',
    'isotype':'IgG',
    'assay_precision':'Quantitative',
    'lab_software':'APRIL TO ASK ABOUT THIS',
    'lab_software_version':'APRIL TO ASK ABOUT THIS',
}

outputs = sdmc_tools.processing_minus_ldms(
    data,
    metadata_dict=md,
    input_data_path=path2,
    cols_to_lower=True,
)

reorder = [
    'network',
    'protocol',
    'specrole',
    'sample_id_lab',
    'ptid',
    'visitno',
    'site_id',
    'drawdt',
    'spectype',
    'upload_lab_id',
    'assay_lab_name',
    'assay_type',
    'assay_subtype',
    'isotype',
    'analyte',
    'assay_precision',
    'instrument',
    'lab_software',
    'lab_software_version',
    'result',
    'result_units',
    'lloq',
    'uloq_analytical',
    'uloq_dynamic',
    'comments_lab',
    'sdmc_processing_datetime',
    'sdmc_data_receipt_datetime',
    'input_file_name',
]
assert set(reorder).symmetric_difference(outputs.columns) == set()

outputs = outputs[reorder]

outputs.to_csv(
    "/networks/vtn/lab/SDMC_labscience/assays/Binding_Ab_via_MSD/PPD_ThermoFisher/example_data/example_processed_data/RRPV-24-02_PPDVAC167_UBT_20260921_v1_DRAFT_PROCESSED_2026-10-05.txt",
    sep="\t",
    index=False
)