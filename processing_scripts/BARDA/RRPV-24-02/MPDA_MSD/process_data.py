import pandas as pd
import numpy as np
import os
from datetime import date

import sdmc_tools.process as sdmc_tools
import sdmc_tools.access_ldms as access_ldms

path2 = '/networks/vtn/lab/SDMC_labscience/assays/Binding_Ab_via_MSD/PPD_ThermoFisher/example_data/RRPV-24-02_PPDVAC167_UBT_20260921_v1.csv'
ubt = pd.read_csv(path2)


## checks against DTA

ubt.PROTOCOL.unique()

# barcode length does not match length in 'data example' from DTA, but i wouldnt have necessarily expected this length to match
ubt.BARCODE.unique()

ubt.SUBJECT.unique()

# subject should match site + 4-digit number
assert np.all(ubt.SUBJECT.str[:3]==ubt.SITE.astype(str))

# visit codes match DTA
assert np.all(ubt.VISIT.isin(['V1D1', 'V2D31', 'V3D91', 'V4D181', 'V5D366', 'ET', 'Unscheduled [X]']))

ubt.COLLECT_DATE.unique()

# check datetime format
parsed = pd.to_datetime(ubt.COLLECT_DATE, format="%Y-%m-%d", errors="coerce")
assert parsed.notna().all()

# check limits of quantitation

limits_from_dta = pd.DataFrame({
    'ANALYTE':[
        'VAC167_ECL_SER_S',
        'VAC167_ECL_SER_LB.1',
        'VAC167_ECL_SER_N',
        'VAC167_ECL_SER_BA.2.86',
        'VAC167_ECL_SER_XBB.1.5v2',
        'VAC167_ECL_SER_KP.3',
        'VAC167_ECL_SER_JN.1v2',
        'VAC167_ECL_SER_KP.2',
        'VAC167_ECL_SER_B.1.617.2',
        'VAC167_ECL_SER_BA.5',
    ],
    'LLOQ':[
        420,
        465,
        725,
        690,
        875,
        570,
        665,
        615,
        740,
        525,
    ],
    'Analytical ULOQ':[
        692599,
        290681,
        1006443,
        371214,
        366762,
        362655,
        370831,
        382728,
        486991,
        305169,
    ],
    'Dynamic ULOQ':[
        63500000,
        11250000,
        44250000,
        16800000,
        21350000,
        13800000,
        16150000,
        14900000,
        45050000,
        31750000,
    ]
})
limits_check = ubt[['ANALYTE','LLOQ','Analytical ULOQ', 'Dynamic ULOQ']].drop_duplicates().merge(
    limits_from_dta,
    on='ANALYTE',
    suffixes=('_data', '_dta')
)
assert len(limits_check) == len(limits_from_dta)
assert (limits_check.LLOQ_data!=limits_check.LLOQ_dta).sum() == 0
assert (limits_check['Analytical ULOQ_data']!=limits_check['Analytical ULOQ_dta']).sum() == 0
assert (limits_check['Dynamic ULOQ_data']!=limits_check['Dynamic ULOQ_dta']).sum() == 0

ubt.RESULT.unique()

ubt.loc[(ubt.RESULT.notna()) & (ubt.RESULT.str.contains(">")), ['RESULT','LLOQ','Analytical ULOQ', 'Dynamic ULOQ']].drop_duplicates()

ubt.loc[(ubt.RESULT.notna()) & (ubt.RESULT.str.contains("<")), ['RESULT','LLOQ','Analytical ULOQ', 'Dynamic ULOQ']].drop_duplicates()

#The DTA does not say that 'nan' is a valid result value, but 'comments' says insufficient sample quantity so i have no questions
ubt.loc[(ubt.RESULT.isna()),]

ubt.RESULT_UNIT.unique()

ubt.loc[ubt.RESULT_UNIT.isna()]

ubt.COMMENTS.unique()


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