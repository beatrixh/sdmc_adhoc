## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 11/19/2025
# Purpose: Accelevir IPDA
## ---------------------------------------------------------------------------##
import pandas as pd
import numpy as np
import os
import datetime as datetime

import sdmc_tools.process as sdmc_tools
import sdmc_tools.access_ldms as access_ldms

def main():
    ipda_path = '/trials/vaccine/p805/s001/qdata/LabData/IPDA_pass-through/20251113 - ACDX IPDAv2 Data Report - HVTN805 cumulative.xlsx'
    ipda = pd.read_excel(ipda_path)

    ipda.columns = [i.lower().replace(" ","_") for i in ipda.columns]

    guspecs = ipda.global_id.str.split(":", expand=True)
    guspecs.columns = ['guspec','guspec2','guspec3']
    ipda = pd.concat([guspecs, ipda], axis=1)

    col_rename = {
        'pbmc_count_at_thaw\xa0':'cell_count_at_thaw',
        'pbmc_viability':'cell_viability',
        'total_diploid_genomes_per_assay_(ipdav2.0)':'total_diploid_genomes',
        'frequency_of_intact_hiv-1_proviruses_per_million_cells_(ipdav2.0)':'result_intact_proviruses',
        'frequency_of_defective_proviruses_per_million_cells_(ipdav2.0)':'result_defective_proviruses',
        'frequency_of_total_hiv-1_proviruses_detected_per_million_cells_(ipdav2.0)':'result_total_proviruses',
        'sample_comments':'comments'
    }
    ipda = ipda.rename(columns=col_rename)

    ipda.sample_quality_issue = ipda.sample_quality_issue.map({"N":False})
    ipda.useit = ipda.useit.map({"Y":True})
    ipda.loc[ipda.signalfail.isna(), 'signalfail'] = False
    ipda.loc[ipda.deviation_of_testing.isna(), 'deviation_of_testing'] = False

    ldms = access_ldms.pull_one_protocol('hvtn', 805)
    assert set(ipda.guspec).union(ipda.guspec2).union(ipda.guspec3).difference(ldms.guspec).difference({None}) == set()
    guspecs_all = list(set(ipda.guspec).union(ipda.guspec2).union(ipda.guspec3).difference({None}))
    ldms = ldms.loc[ldms.guspec.isin(guspecs_all)]

    metadata = pd.read_excel('/networks/vtn/lab/SDMC_labscience/assays/IPDA/Accelevir/SDMC_materials/IPDA_info.xlsx', header=None)

    # we'll populate these ourselves
    ipda = ipda.drop(columns=[
        'assay_lab_name',
        'assay_name',
    ])

    # standard processing
    md = metadata.set_index(0)[1].to_dict()
    md['Specrole'] = 'Sample'

    outputs = sdmc_tools.standard_processing(
        input_data=ipda,
        input_data_path=ipda_path,
        guspec_col='guspec',
        network='hvtn',
        metadata_dict=md,
        ldms=ldms
    )

    outputs.collection_date_time = pd.to_datetime(outputs.collection_date_time).dt.date.astype(str)
    outputs.drawdt = outputs.drawdt.astype(str)

    outputs.ptid = outputs.ptid.astype(int)
    outputs.visitno = outputs.visitno.astype(float)
    outputs.protocol = outputs.protocol.astype(float)

    assert (outputs.drawdt!=outputs.collection_date_time).sum() == 0
    assert (outputs.ptid!=outputs.patient_id).sum() == 0
    assert (outputs.visitno!=outputs.visit_id).sum() == 0
    assert (outputs.collection_date_time!=outputs.drawdt).sum() == 0

    assert len(outputs[['protocol','study_id']].drop_duplicates()) == 1

    outputs = outputs.drop(columns=[
        'study_id',
        'patient_id',
        'visit_id',
        'collection_date_time',
        'spec_type',
        'global_id',
    ])
    outputs = outputs.rename(columns={'guspec':'guspec1'})

    reorder = [
        'network',
        'protocol',
        'specrole',
        'guspec1',
        'guspec2',
        'guspec3',
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
        'assay_subtype',
        'assay_precision',
        'instrument',
        'assay_date',
        'cell_count_at_thaw',
        'cell_viability',
        'cd4+_t_cell_count',
        'cd4+_t_cell_viability',
        'dna_shearing_index',
        'total_diploid_genomes',
        'result_intact_proviruses',
        'result_defective_proviruses',
        'result_total_proviruses',
        'result_units',
        'comments',
        'sample_quality_issue',
        'deviation_of_testing',
        'useit',
        'signalfail',
        'sdmc_processing_datetime',
        'sdmc_data_receipt_datetime',
        'input_file_name',
    ]
    assert set(outputs.columns).symmetric_difference(reorder) == set()
    outputs = outputs[reorder]

    savedir = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN805_HPTN093/assays/IPDA/misc_files/data_processing/'
    today = datetime.date.today().isoformat()
    outputs.to_csv(
        savedir + f"HVTN805_Accelevir_IPDA_Processed_{today}.txt",
        sep="\t",
        index=False
    )

    summary = pd.pivot_table(
        outputs,
        index=['ptid'],
        columns='visitno',
        aggfunc='count',
        values='result_total_proviruses',
        fill_value=0
    )

    summary.to_excel(
        savedir + "HVTN805_Accelevir_IPDA_summary.xlsx"
    )

    manifest = pd.read_excel(
        '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN805_HPTN093/assays/IPDA/misc_files/MANIFEST_SHIP_REQ_0352.xlsx',
    )

    assert set(manifest.GLOBAL_ID).symmetric_difference(guspecs_all)==set()

if __name__=="__main__":
    main()