## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 01.29.2026
# Purpose: MTB Rif: Draft -- reshape example data into standard format
## ---------------------------------------------------------------------------##

import pandas as pd
import numpy as np
import sdmc_tools.access_ldms as access_ldms
import sdmc_tools.process as sdmc
import os
import datetime

from io import StringIO

example_path = "/networks/vtn/lab/SDMC_labscience/assays/MTB-Rif/example_data/CTAC0127_2_2025.12.08_16.02.28.csv"
with open(example_path,"r", encoding='utf-16') as f:
    data = f.read()

data = data.split('\n')

result_idx = [i for i in range(len(data)) if data[i]=='RESULT TABLE']

# separate data out into all the rows that correspond to each sample
result_blocks = {}
for i in result_idx:
    sample_id = data[int(i) + 1].split(",")[1]
    j = i
    while data[j]!="Test stopped by decision command [1].":
        j += 1
        result_blocks[sample_id] = data[i:j]

def make_dict(d):
    out = {}
    for i in d:
        i = i.split(",")
        if len(i)>1:
            out[i[0]] = [i[1]]
    return out
    

# for each sample, pull out cartridge/analyte/result sub-data
result_sub_blocks = {}
for sample, block in result_blocks.items():
    result_sub_blocks[sample] = {}

    j = 0
    # result table looks like a dict
    while block[j]!="":
        j += 1
        result_table = block[:j]
    result_table = pd.DataFrame(make_dict(result_table))
    result_sub_blocks[sample]['result_table'] = result_table

    # cartidge index looks like a dict
    while block[j]!='Cartridge Index,1':
        j +=1
    k = j
    while block[k]!="":
        k += 1
        cartridge_index = block[j:k]
    cartridge_index = pd.DataFrame(make_dict(cartridge_index))
    result_sub_blocks[sample]['cartridge_index'] = cartridge_index

    # analyte result looks like a csv
    while block[k]!='Analyte Result':
        k +=1
    l = k
    while block[l]!="":
        l += 1
        analyte_result = block[k+1:l]
    analyte_result = pd.read_csv(StringIO("\n".join(analyte_result)))
    result_sub_blocks[sample]['analyte_result'] = analyte_result

# concat everything into one table
df = pd.DataFrame()
for k, d in result_sub_blocks.items():
    # print(k)
    sub = pd.merge(
        d['result_table'],
        d['analyte_result'],
        how='cross'
    ).merge(
        d['cartridge_index'],
        how='cross'
    )
    df = pd.concat([df, sub])

savepath="/networks/vtn/lab/SDMC_labscience/assays/MTB-Rif/SDMC_materials/EXAMPLE_PROCESSED_CTAC0127_2_2025.12.08_16.02.28.csv"

df.to_csv(savepath, sep="\t", index=False)