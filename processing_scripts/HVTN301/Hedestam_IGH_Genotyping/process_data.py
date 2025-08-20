## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 08/20/2025
# Purpose: Process IGH Genotyping data from Hedestam Lab (Karolinka)
## ---------------------------------------------------------------------------##
import pandas as pd
import numpy as np
import os
import datetime as datetime

import sdmc_tools.process as sdmc_tools
import sdmc_tools.access_ldms as access_ldms

def main():
    # read in inputs data
    input_data_path = "/trials/vaccine/p301/s001/qdata/LabData/IGHV_genotyping_pass-through/HVTN301_IGH_GENOTYPES.xlsx"

    df = pd.read_excel(
        input_data_path,
        sheet_name=None
    )

    # check columns are identical
    assert len(set(df['IGHJ'].columns).symmetric_difference(df['IGHD'].columns)) == 0
    assert len(set(df['IGHJ'].columns).symmetric_difference(df['IGHV'].columns)) == 0

    # cast to long
    data = []
    for k in df.keys():
        d = df[k]
        d = pd.melt(
            d,
            value_vars=list(d.columns),
            var_name='ptid',
            value_name='allele'
        )
        d['region'] = k
        d['gene'] = d.allele.str.split('*', expand=True)[0]
        data += [d]

    data = pd.concat(data)
    data = data.loc[data.allele.notna()]

    # cast to wide on allele 1/2/3/4, leave 'na' for ptid/genes with < 4 alleles
    data['allele_id'] = 1
    data['allele_id'] = 'allele_' + data.groupby(['ptid','gene']).allele_id.transform('cumsum').astype(str)

    data = pd.pivot(data, index=['ptid','region','gene'], columns='allele_id', values='allele')
    data = data.reset_index()
    data.columns.name = ''

    # check against manifest
    mfst = pd.read_csv('/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN301/assays/genotyping/misc_files/shipping_manifest.txt', sep="\t")
    assert len(set(data.ptid).symmetric_difference(mfst.PID))==0

    # make sure there are exactly 53 of each of these
    mfst[['PID','GLOBAL_ID']].nunique()

    # merge guspec onto data using ptid
    mfst = mfst.rename(columns={'GLOBAL_ID':'guspec', 'PID':'ptid'})
    data = data.merge(mfst[['guspec','ptid']], on='ptid', how='left')

    # oull ldms
    ldms = access_ldms.pull_one_protocol('hvtn', 301)
    ldms = ldms.loc[ldms.guspec.isin(data.guspec)]

    # standard processing
    md = {
        'network':'HVTN',
        'upload_lab_id': 'K6',
        'assay_lab_name':'Hedestam Lab (Karolinska)',
        'assay_type': 'IGH Genotyping',
        'assay_subtype':'2x300bp',
        'instrument':'Illumina MiSeq',
        'lab_software':'IgDiscover'
    }

    outputs = sdmc_tools.standard_processing(
        input_data=data,
        input_data_path=input_data_path,
        guspec_col='guspec',
        network='hvtn',
        metadata_dict=md,
        ldms=ldms
    )

    # make sure ptid-guspec combos from manifest match those in LDMS
    assert (outputs.ptid_x.astype(int)!=outputs.ptid_y.astype(int)).sum()==0
    outputs = outputs.drop(columns='ptid_y').rename(columns={'ptid_x':'ptid'})

    # reorder columns
    reorder = [
        'network',
        'protocol',
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
        'assay_subtype',
        'instrument',
        'lab_software',
        'region',
        'gene',
        'allele_1',
        'allele_2',
        'allele_3',
        'allele_4',
        'sdmc_processing_datetime',
        'sdmc_data_receipt_datetime',
        'input_file_name'
    ]

    assert len(set(reorder).symmetric_difference(outputs.columns)) == 0
    outputs = outputs[reorder]

    # save to .txt
    savedir = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN301/assays/genotyping/misc_files/data_processing/'
    today = datetime.date.today().isoformat()
    outputs.to_csv(savedir + f"HVTN301_IGH_Genotyping_processed_{today}.txt", sep="\t", index=False)

    # pivot summaries
    pivot_summary = pd.melt(
        outputs,
        id_vars=['guspec','ptid','visitno','region','gene'],
        value_vars=['allele_1','allele_2','allele_3','allele_4'],
        value_name='allele'
    )

    pivot_summary = pivot_summary.loc[pivot_summary.allele.notna()]
    summary1 = pd.pivot_table(
        pivot_summary,
        index=['guspec','ptid','visitno'],
        columns=['region'],
        aggfunc='count'
    )[['allele']]

    summary2 = pd.pivot_table(
        pivot_summary,
        index=['guspec','ptid'],
        columns=['visitno'],
        aggfunc='count',
        fill_value=0
    )[['allele']]

    summary1.to_excel(savedir + "HVTN301_IGH_Genotyping_region_summary.xlsx")
    summary2.to_excel(savedir + "HVTN301_IGH_Genotyping_visitno_summary_.xlsx")

if __name__=="__main__":
    main()
