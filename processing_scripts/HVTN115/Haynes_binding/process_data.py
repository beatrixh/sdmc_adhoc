## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 11/12/2025
# Purpose: Binding
## ---------------------------------------------------------------------------##
import pandas as pd
import numpy as np
import os
import datetime as datetime

import sdmc_tools.process as sdmc_tools
import sdmc_tools.access_ldms as access_ldms

# read in binding
# read in controls
# merge on controls
# subset to 115/135
    # standard formatting
    # save etc

def main():
    # read in data ------------------------------------------------------------##
    binding_path = '/trials/vaccine/p115/s001/qdata/binding_blocking_pass-through/250930 HVTN115 and 135 post-5 binding summary.xlsx'
    binding = pd.read_excel(binding_path, sheet_name='summary', skiprows=5)

    # these should uniquely identify a row
    assert binding[['Subject ID','Ag letter']].drop_duplicates().shape[0] == len(binding)

    # cast to long, formatting
    binding = binding.melt(
        id_vars=['Run Date', 'Plate ', 'sample position', 'Ag letter', 'Ag Lot',
           'Original ID', 'Date Drawn', 'Infant-Mother Pair Match # (HVTN135)',
           'Study ID', 'Group', 'Tx Group',
           'Part (HVTN115) or subject type (HVTN135)', 'Pub ID', 'Subject ID',
           'Visit', 'Day', 'Ag','Log AUC', 'EC50', 'End Titer'],
        value_vars=['_30', '_90', '_270', '_810', '_2430', '_7290',
           '_21870', '_65610', '_196830', '_590490', '_1771470', '_5314410'],
        value_name='result',
        var_name='dilution'
    )

    binding.dilution = binding.dilution.str[1:].astype(int)
    binding['specrole'] = 'Sample'

    binding = binding.rename(columns={
        'Plate ':'plate_number',
        'Original ID':'guspec',
    })
    binding.columns = [i.lower().replace(" ","_") for i in binding.columns]

    # read in controls --------------------------------------------------------##
    controls = pd.read_excel(binding_path, sheet_name='controls', skiprows=3)

    controls = controls.melt(
        id_vars=['Ag letter', 'Plate #', 'Ag lot', 'Ag', 'Control','*Log AUC'],
        value_vars=['_100',
           '_33.333', '_11.1111', '_3.7037', '_1.2346', '_0.4115', '_0.1372',
           '_0.0457', '_0.0152', '_0.0051', '_0.001697', '_0.000564',],
        value_name='result',
        var_name='control_concentration'
    )
    controls = controls.rename(columns={'Control':'control_name'})
    controls.control_concentration = controls.control_concentration.str[1:].astype(float)
    controls['control_concentration_units'] = 'ug/ml'
    controls['specrole'] = 'Standard'

    controls = pd.concat([
        pd.DataFrame(data={
        'specrole':['Background'],
        'result':[0.044],}),
        controls
    ])

    controls.columns = [i.lower().replace(" ","_") for i in controls.columns]

    controls = controls.rename(columns={
        'plate_#':'plate_number',
        '*log_auc':'log_auc',
    })

    data = pd.concat([binding, controls])

    # standard formatting -----------------------------------------------------##
    ldms = access_ldms.pull_multiple_protocols('hvtn', [115, 135])
    ldms = ldms.loc[ldms.guspec.isin(data.guspec.tolist())]

    # standard processing
    md = {
        'upload_lab_id': 'BH',
        'assay_lab_name':'Haynes Lab (Duke)',
        'assay_type': 'ELISA',
        'assay_subtype':'Binding',
        'instrument':'SpectraMax',
        'lab_software':'Softmax',
        'lab_software_version':'Softmax 5.3',
        'result_units':'Absorbance (450 nm)',
        'assay_precision':'Quantitative',
    }

    outputs = sdmc_tools.standard_processing(
        input_data=data,
        input_data_path=binding_path,
        guspec_col='guspec',
        network='hvtn',
        metadata_dict=md,
        ldms=ldms
    )

    assert len(outputs.loc[(outputs.ptid.notna()) & (outputs.ptid.astype(float)!=outputs.subject_id.astype(float)),['ptid','subject_id']]) == 0

    outputs.loc[outputs.guspec.isna(),['specrole']].value_counts()

    # misc checks -------------------------------------------------------------##
    # visual check -- for non-na, do these match? yes
    outputs[['study_id','protocol']].drop_duplicates()

    # for rows w/out guspec, map protocol from lab submitted data
    outputs = outputs.drop(columns=['study_id'])

    checkrows = outputs.loc[(outputs.specrole=='Sample')]

    assert (checkrows.ptid.astype(float)!=checkrows.subject_id.astype(float)).sum() == 0
    assert (checkrows.visitno.astype(float)!=checkrows.visit.astype(float)).sum() == 0
    assert len(checkrows.loc[(checkrows.date_drawn!=checkrows.drawdt) & (checkrows.date_drawn.notna()),['date_drawn','drawdt']]) == 0
    assert len(checkrows.loc[checkrows.ptid.astype(float)!=checkrows.subject_id.astype(float), ['ptid','subject_id']]) == 0

    # reformatting ------------------------------------------------------------##
    outputs = outputs.drop(columns=['date_drawn','subject_id','visit'])

    outputs = outputs.rename(columns={
        'ag':'analyte',
        'ag_lot':'analyte_lot'
    })

    outputs = outputs.drop(columns=[
        'pub_id',
        'day',
        'infant-mother_pair_match_#_(hvtn135)',
        'part_(hvtn115)_or_subject_type_(hvtn135)',
        'group',
        'tx_group',
        'ag_letter',
    ])

    outputs.loc[outputs.specrole=='Sample','network'] = 'HVTN'

    reorder = [
        'network',
        'protocol',
        'specrole',
        'guspec',
        'ptid',
        'visitno',
        'drawdt',
        'run_date',
        'spectype',
        'spec_primary',
        'spec_additive',
        'spec_derivative',
        'upload_lab_id',
        'assay_lab_name',
        'assay_type',
        'assay_subtype',
        'assay_precision',
        'instrument',
        'lab_software',
        'lab_software_version',
        'analyte',
        'analyte_lot',
        'control_name',
        'control_concentration',
        'control_concentration_units',
        'dilution',
        'result',
        'result_units',
        'ec50',
        'end_titer',
        'log_auc',
        'plate_number',
        'sample_position',
        'sdmc_processing_datetime',
        'sdmc_data_receipt_datetime',
        'input_file_name',
    ]

    assert set(outputs.columns).symmetric_difference(reorder) == set()
    outputs = outputs[reorder]

    # additional checks -------------------------------------------------------##
    ## check we still have expected number of ag-ptids

    assert len(outputs.loc[outputs.specrole=='Sample',['ptid','analyte']].drop_duplicates()) == 686

    # # visual check, ptid/analyte uniquely identifies 686 unique combos
    # outputs.loc[outputs.specrole=='Sample'].groupby(['ptid','analyte']).count().result.value_counts()

    # visual check, either sample with protocol == 115/135, or background/standard
    outputs[['specrole','protocol']].drop_duplicates()

    # visual check, missing guspecs only for non-samples
    outputs.loc[outputs.guspec.isna(),['protocol','specrole']].drop_duplicates()

    # visual check
    outputs.loc[outputs.result.isna()].control_name.value_counts()

    # sara confirmed we should drop these
    outputs = outputs.loc[outputs.result.notna()]

    # save data and pivot summaries -------------------------------------------##
    outputs.protocol = outputs.protocol.astype(float)

    outputs115 = outputs.loc[(outputs.protocol==115.) | (outputs.specrole.isin(['Standard','Background']))]
    outputs135 = outputs.loc[(outputs.protocol==135.) | (outputs.specrole.isin(['Standard','Background']))]

    today = datetime.date.today().isoformat()

    # 115
    savedir115 = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN115/assays/binding_ELISA/misc_files/data_processing/'
    outputs115.to_csv(savedir115 + f'HVTN115_ELISA_binding_processed_{today}.txt', sep="\t", index=False)
    summary115 = pd.pivot_table(
        outputs115,
        index=['specrole','ptid','visitno'],
        columns=['analyte'],
        values='result',
        aggfunc='count',
        dropna=False
    )
    summary115.dropna(how='all').to_excel(savedir115 + "HVTN115_ELISA_bnab_binding_summary.xlsx")

    # 135
    savedir135 = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN135/assays/binding_ELISA/misc_files/data_processing/'
    outputs135.to_csv(savedir135 + f'HVTN135_ELISA_binding_processed_{today}.txt', sep="\t", index=False)
    summary135 = pd.pivot_table(
        outputs135,
        index=['specrole','ptid','visitno'],
        columns=['analyte'],
        values='result',
        aggfunc='count',
        dropna=False,
    )
    summary135.dropna(how='all').to_excel(savedir135 + "HVTN135_ELISA_bnab_binding_summary.xlsx")

    ## final checks -----------------------------------------------------------##

    binding_path = '/trials/vaccine/p115/s001/qdata/binding_blocking_pass-through/250930 HVTN115 and 135 post-5 binding summary.xlsx'
    raw = pd.read_excel(binding_path, sheet_name='summary', skiprows=5)

    set(raw.loc[raw['Study ID']=='HVTN 135']['Subject ID'].unique()).symmetric_difference(outputs135.ptid.dropna())

    set(raw.loc[raw['Study ID']=='HVTN115']['Subject ID'].unique()).symmetric_difference(outputs115.ptid.dropna())

    # completeness check
    manifest1 = pd.read_csv(
        '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN115/assays/binding_ELISA/misc_files/manifests/481-373-0000043537.txt',
        sep='\t'
    )
    manifest2 = pd.read_csv(
        '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN115/assays/binding_ELISA/misc_files/manifests/512--2025000186.txt',
        sep='\t'
    )
    manifest = pd.concat([manifest1, manifest2])

    # partial match :-/
    set(manifest.loc[manifest.PROTOCOL==135.].GLOBAL_ID).symmetric_difference(outputs135.guspec)

    # full match
    set(manifest.loc[manifest.PROTOCOL==115.].GLOBAL_ID).symmetric_difference(outputs115.guspec)

if __name__=="__main__":
    main()