## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Save list of novel alleles not in novel allele file
## ---------------------------------------------------------------------------##
import pandas as pd
import numpy as np
import os
import datetime as datetime

import sdmc_tools.process as sdmc_tools
import sdmc_tools.access_ldms as access_ldms

# save novel allele sequences to fasta file
datadir = '/trials/vaccine/p144/s001/qdata/LabData/genotyping_pass-through/'
novels_h = pd.read_excel(datadir + 'GENOTYPES_IGH_FILT.xlsx', sheet_name='Novel_alleles').rename(columns={'IGHV':0})
novels_k = pd.read_excel(datadir + 'GENOTYPES_IGK_FILT.xlsx', sheet_name='IGKV_novel_alleles', header=None)
novels_l = pd.read_excel(datadir + 'GENOTYPES_IGL_FILT.xlsx', sheet_name='IGLV_Novel_Alleles', header=None)

novels = pd.concat([novels_h, novels_k, novels_l])
novel_in_extra = novels.loc[novels[0].str.contains(">")][0].str[1:].unique()

# read in 301 data
input_data_path = "/trials/vaccine/p301/s001/qdata/LabData/IGH_genotyping_pass-through/HVTN301_IGH_GENOTYPES.xlsx"

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
novel301 = data.loc[data.allele.str.rpartition("_")[2].str[0]=="S"].allele.unique()

# read in 144
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

    # # cast to wide on allele 1/2/3/4, leave 'na' for ptid/genes with < 4 alleles
    # data['allele_id'] = 1
    # data['allele_id'] = 'allele_' + data.groupby(['ptid','gene']).allele_id.transform('cumsum').astype(str)

    # data = pd.pivot(data, index=['ptid','region','gene'], columns='allele_id', values='allele')
    # data = data.reset_index()
    # data.columns.name = ''
    return data

datadir = '/trials/vaccine/p144/s001/qdata/LabData/genotyping_pass-through/'
fnames = [i for i in os.listdir(datadir) if '.xlsx' in i and '$' not in i]
datasets = [pd.read_excel(datadir + i, sheet_name=None) for i in fnames]

data = pd.concat([format_data(d) for d in datasets])
novel144 = data.loc[data.allele.str.rpartition("_")[2].str[0]=="S"].allele.unique()

novel144 = pd.DataFrame({'alleles in 144':list(set(novel144).difference(novel_in_extra))})
novel301 = pd.DataFrame({'alleles in 301':list(set(novel301).difference(novel_in_extra))})

savedir = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN144/assays/genotyping/misc_files/data_processing/'
# Save to Excel with separate tabs
with pd.ExcelWriter(savedir + "additional_novel_alleles.xlsx", engine="xlsxwriter") as writer:
    novel144.to_excel(writer, sheet_name="alleles in 144", index=False)
    novel301.to_excel(writer, sheet_name="alleles in 301", index=False)