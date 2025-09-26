## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 09/26/2025
# Purpose: Process Ig Genotyping data from Hedestam Lab (Karolinka)
## ---------------------------------------------------------------------------##
import pandas as pd
import numpy as np
import os
import datetime as datetime

import sdmc_tools.process as sdmc_tools
import sdmc_tools.access_ldms as access_ldms

def main():
    datadir = '/trials/vaccine/p144/s001/qdata/LabData/genotyping_pass-through/'

    # cast to long
    def format_data(excel_dict):
        data = []
        sheets_to_use = [i for i in excel_dict.keys() if len(i)==4] #should be "IGHV", "IGHJ", etc
        for k in sheets_to_use:
            d = excel_dict[k]
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
        return data

    datadir = '/trials/vaccine/p144/s001/qdata/LabData/genotyping_pass-through/'
    datasets = [pd.read_excel(datadir + i, sheet_name=None) for i in os.listdir(datadir)]

    data = pd.concat([format_data(d) for d in datasets])

    mfst1 = pd.read_csv('/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN144/assays/genotyping/misc_files/512--2024000249.txt', sep="\t")
    mfst2 = pd.read_csv('/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN144/assays/genotyping/misc_files/512--2025000058.txt', sep="\t")


    mfst = pd.concat([mfst1, mfst2])
    # set(mfst.PID).symmetric_difference(data.ptid)
    # mfst.PID.nunique(), len(mfst)

    data['guspec'] = data.ptid.map(mfst.set_index('PID').GLOBAL_ID)

    ldms = access_ldms.pull_one_protocol('hvtn', 144)
    ldms = ldms.loc[ldms.guspec.isin(data.guspec)]

    # set(ldms.guspec).symmetric_difference(data.guspec)

    # standard processing
    md = {
        'network':'HVTN',
        'upload_lab_id': 'K6',
        'assay_lab_name':'Hedestam Lab (Karolinska)',
        'assay_type': 'Ig Genotyping',
        'assay_subtype':'2x300bp',
        'instrument':'Illumina MiSeq',
        'lab_software':'IgDiscover',
        'lab_software_version':'IgDiscover22 v1.0.5.',
    }

    outputs = sdmc_tools.standard_processing(
        input_data=data.drop(columns='ptid'),
        input_data_path=datadir,
        guspec_col='guspec',
        network='hvtn',
        metadata_dict=md,
        ldms=ldms
    )

    outputs['input_file_name'] = 'GENOTYPES_IG'+outputs.region.str[2]+'_FILT.xlsx'

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
        'lab_software_version',
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
    savedir = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN144/assays/genotyping/misc_files/data_processing/'
    today = datetime.date.today().isoformat()
    outputs.to_csv(savedir + f"HVTN144_Ig_Genotyping_processed_{today}.txt", sep="\t", index=False)

    # pivot summaries
    pivot_summary = pd.melt(
        outputs,
        id_vars=['guspec','ptid','visitno','region','gene'],
        value_vars=['allele_1','allele_2','allele_3','allele_4'],
        value_name='allele'
    )

    pivot_summary = pivot_summary.loc[pivot_summary.allele.notna()]

    summary2 = pd.pivot_table(
        pivot_summary,
        index=['guspec','ptid'],
        columns=['visitno'],
        aggfunc='count',
        fill_value=0
    )[['allele']]

    pivot_summary['category'] = pivot_summary.region.str[2]
    summary0 = pd.pivot_table(
        pivot_summary,
        index=['guspec','ptid','visitno'],
        columns=['category','region'],
        aggfunc='count',
        fill_value=0
    )['allele']



    summary0.to_excel(savedir + "HVTN144_Ig_Genotyping_region_summary.xlsx")
    summary2.to_excel(savedir + "HVTN144_Ig_Genotyping_visitno_summary_.xlsx")

if __name__=="__main__":
    main()
