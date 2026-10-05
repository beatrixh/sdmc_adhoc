## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 2026-10-05
# Purpose: Process MOCK monogram nab data for BARDA. following sara's spec
## ---------------------------------------------------------------------------##

import pandas as pd
import numpy as np
import os
from datetime import date

import re
from datetime import datetime

import sdmc_tools.process as sdmc_tools
import sdmc_tools.access_ldms as access_ldms

path = '/networks/vtn/lab/SDMC_labscience/assays/NAb/Monogram/example_data/RRPV-24-02_RETAIL_Monogram_YYYYMMDD_NAb_D614G_KP2_updated.csv'
df = pd.read_csv(path)

## TESTS ---------------------------------------------------------------------##
def invalid_rows(df: pd.DataFrame, column: str, is_valid, dropna) -> pd.DataFrame:
    """Return non-empty values in `column` that fail `is_valid`."""
    values = df[column]
    if dropna:
        values = df[column].fillna("")
    mask = (values != "") & ~values.map(is_valid)
    return df.loc[mask, [column]]


def assert_column_exists(df: pd.DataFrame, column: str) -> None:
    assert column in df.columns, f"Required column is missing: {column}"


def assert_no_invalid_values(
    df: pd.DataFrame,
    column: str,
    is_valid,
    rule: str,
    dropna = False,
) -> None:
    assert_column_exists(df, column)
    bad = invalid_rows(df, column, is_valid, dropna)
    assert bad.empty, (
        f"{column} must {rule}. "
        f"Invalid rows: {bad.index.tolist()}; "
        f"Invalid values: {bad[column].tolist()}"
    )

def test_protocol_code(df):
    assert_no_invalid_values(
        df,
        "PROTOCOL_CODE",
        lambda value: value == "RRPV-24-02",
        'equal "RRPV-24-02"',
    )
test_protocol_code(df)

def test_monogram_accession_format(df):
    assert_no_invalid_values(
        df,
        "MONOGRAM_ACCESSION",
        lambda value: bool(re.fullmatch(r"\d{2}-\d{6}", value)),
        'match "XX-XXXXXX" (for example, "26-123456")',
    )

test_monogram_accession_format(df)

def test_sample_id_format(df):
    assert_no_invalid_values(
        df,
        "SAMPLE_ID",
        lambda value: bool(re.fullmatch(r"\d{10}-\d{2}", value)),
        'match "XXXXXXXXXX-XX" (10 digits, hyphen, 2 digits)',
    )
test_sample_id_format(df)

def test_patient_number_format(df):
    assert_no_invalid_values(
        df,
        "PATIENT_NUMBER",
        lambda value: bool(re.fullmatch(r"\d{3}-\d{4}", value)),
        'match "###-####"',
    )
test_patient_number_format(df)

def test_collection_date_format(df):
    def is_valid_date(value: str) -> bool:
        try:
            datetime.strptime(value, "%d-%b-%Y")
            return True
        except ValueError:
            return False

    assert_no_invalid_values(
        df,
        "COLLECTION_DATE",
        is_valid_date,
        'be a real date in "dd-Mon-yyyy" format (for example, "04-Sep-2026")',
    )
test_collection_date_format(df)

def test_collection_time_format(df):
    def is_valid_time(value: str) -> bool:
        try:
            datetime.strptime(value, "%H:%M")
            return True
        except ValueError:
            return False

    assert_no_invalid_values(
        df,
        "COLLECTION_TIME",
        is_valid_time,
        'be a valid 24-hour time in "hh:mm" format (for example, "13:45")',
    )
test_collection_time_format(df)

def test_pvc(df):
    pvc_allowable_values = ['V1D1', 'V2D31', 'V3D91', 'V4D181', 'V5D366', 'ET']
    assert_no_invalid_values(
        df,
        "PVC",
        lambda value: value in pvc_allowable_values,
        "be included in pvc_allowable_values",
    )
test_pvc(df)

def test_sars_cov_2_variant(df):
    allowed = ["SARS-COV-2 D614G", "SARS-COV-2 KP2"]

    # Replace this with the exact column name in your dataframe if needed.
    column = "TEST_METHOD"

    assert_no_invalid_values(
        df,
        column,
        lambda value: value in allowed,
        f"be one of {sorted(allowed)}",
    )
test_sars_cov_2_variant(df)

df.TEST_METHOD.unique()

# # STILL HAS LEADING SPACES
# df.TEST_METHOD = df.TEST_METHOD.str.strip()
# test_sars_cov_2_variant(df)

# df = df.rename(columns={'TEST METHOD':'TEST_METHOD'})

# df.TEST_METHOD = ["SARS-COV-2 D614G", "SARS-COV-2 KP.2","SARS-COV-2 D614G", "SARS-COV-2 KP.2","SARS-COV-2 D614G", "SARS-COV-2 KP.2"] ## didnt pass this one

def test_titer_result(df):
    titer_result_allowable_values = ['<40', 'ZNG40', 'NRR']
    def is_numeric_or_allowable(value: str) -> bool:
        if value in titer_result_allowable_values:
            return True
        try:
            float(value)
            return True
        except ValueError:
            return False

    assert_no_invalid_values(
        df,
        "Titer Result",
        is_numeric_or_allowable,
        "be numeric or be included in titer_result_allowable_values",
    )
test_titer_result(df)

def test_unit_of_results(df):
    assert_no_invalid_values(
        df,
        "Unit of Results",
        lambda value: value == "1/Dilution",
        'equal "1/Dilution"',
    )
test_unit_of_results(df)

def test_reason_not_done_or_comment_for_result(df):
    reason_not_done_allowable_values = ['ZPABG', 'QNS']
    assert_no_invalid_values(
        df,
        "Reason Not Done or Comment for Result",
        lambda value: value in reason_not_done_allowable_values,
        "be included in reason_not_done_allowable_values",
        dropna=True
    )
test_reason_not_done_or_comment_for_result(df)

def test_specimen_type(df):
    assert_no_invalid_values(
        df,
        "SPECIMEN_TYPE",
        lambda value: value == "Serum",
        'equal "Serum"',
    )
test_specimen_type(df)

df.TEST_METHOD.unique()

## STANDARD PROCESSING

# merge limits on from DTP # this is the version we want to use
# limits = pd.DataFrame({
#     'TEST_METHOD':['SARS-COV-2 D614G','SARS-COV-2 KP2'],
#     'lloq':[52, 48],
#     'uloq':[93447, 28437]
# })
# df = df.merge(limits, on = 'TEST_METHOD')

# merge limits on from DTP # version we are using with leading spaces
limits = pd.DataFrame({
    'TEST_METHOD':[' SARS-COV-2 D614G',' SARS-COV-2 KP2'],
    'lloq':[52, 48],
    'uloq':[93447, 28437]
})
df = df.merge(limits, on = 'TEST_METHOD')


# add metadata
md = {
    'network':'BARDA',
    'specrole':'Sample',
    'upload_lab_id':'L2',
    'assay_lab_name':'Monogram Biosciences',
    'assay_type':'Neutralizing Antibody',
    'assay_subtype':'PhenoSense SARS CoV-2',
    'assay_details':'HEK293 target cells',
    'assay_precision':'Quantitative',
    'cutoff':50,
}



outputs = sdmc_tools.processing_minus_ldms(
    df,
    metadata_dict=md,
    input_data_path=path,
    cols_to_lower=True,
)

# renaming, following sara's spec
rename = {
    'titer_result':'result_titer',
    'unit_of_results':'result_units',
    'reason_not_done_or_comment_for_result':'lab_comments',
    'collection_date':'drawdt',
    'monogram_accession':'sample_id_lab',
    'specimen_type':'spectype',
    'protocol_code':'protocol',
    'pvc':'visitno',
    'patient_number':'ptid',
}
outputs = outputs.rename(columns=rename)
outputs['pseudovirus'] = outputs.test_method.str[11:]

# expecting to get this column from their manifest in the real data
outputs['visit_description'] = 'MISSING_FROM_MOCK_DATA'

# reorder according to sara's spec
reorder = [
    'network',
    'protocol',
    'sample_id',
    'sample_id_lab',
    'specrole',
    'visitno',
    'visit_description',
    'ptid',
    'drawdt',
    'collection_time',
    'spectype',
    'upload_lab_id',
    'assay_lab_name',
    'test_method',
    'assay_type',
    'assay_subtype',
    'assay_details',
    'assay_precision',
    'pseudovirus',
    'result_titer',
    'result_units',
    'lloq',
    'uloq',
    'cutoff',
    'lab_comments',
    'input_file_name',
    'sdmc_data_receipt_datetime',
    'sdmc_processing_datetime',
]


# make sure these match
assert set(reorder).symmetric_difference(outputs.columns) == set()

outputs = outputs[reorder]


# save to .txt --------------------------------------------------------------------##

today = date.today().isoformat()
fname = f"RRPV-24-02_RETAIL_Monogram_YYYYMMDD_NAb_D614G_KP2_updated_processed_by_sdmc_{today}.txt"

savedir1 = f"/networks/vtn/lab/SDMC_labscience/assays/NAb/Monogram/example_data/"
outputs.to_csv(savedir1 + fname, index=False)

savedir2 = f"/networks/vtn/lab/SDMC_labscience/studies/BARDA/RRPV-24-02/assays/nAb/misc_files/data_processing/"
outputs.to_csv(savedir2 + fname, index=False)
