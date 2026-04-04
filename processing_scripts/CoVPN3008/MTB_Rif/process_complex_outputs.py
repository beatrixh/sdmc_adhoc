import pandas as pd
import numpy as np
import sdmc_tools.access_ldms as access_ldms
import sdmc_tools.process as sdmc
import os
import datetime

from io import StringIO

def process_complex_submission(data_path):
    with open(data_path,"r", encoding='utf-16') as f:
        data = f.read()
    
    data = data.split('\n')
    data = [i.strip("\t").strip('"').replace("\t",",") for i in data] # had to add in random things here for handling differences between submissions
    
    result_idx = [i for i in range(len(data)) if 'RESULT TABLE' in data[i]]
    
    # separate data out into all the rows that correspond to each sample
    result_blocks = {}
    for i in result_idx:
        sample_id = data[int(i) + 1].split(",")[1]
        j = i
        while "Test stopped by decision command [1]." not in data[j]:
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
        while 'Cartridge Index,1' not in block[j]:
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

    df = df.rename(columns={'Test Result':'test_result'})
    
    result_detail = df.test_result.str.strip("|").str.replace("||","|").str.split("|", expand=True).rename(columns={0:'result_mtb',1:'result_rif'})
    df = pd.concat([
        df.drop(columns='test_result'),
        result_detail
    ], axis=1)
    
    return df

ldms = access_ldms.pull_one_protocol('covpn', 3008)

md = {
    'network': 'CoVPN',
    'protocol': 3008,
    'specrole':'Sample',
    'upload_lab_id': 'CH',
    'assay_lab_name': 'CHIL',
    'assay_type': 'RT-PCR', # ? Xpert MTB/RIF Ultra?
    'assay_subtype': 'Xpert MTB/RIF Ultra',
    'assay_details':'Assay modified for tongue swab specimens by diluting SR Buffer',
    'buffer_diluent':'Phosphate Buffered Saline (PBS)',
    'buffer_dilution':'66%',
    'instrument': 'P049 Xpert MTB/RIF Ultra', # ?
    'lab_software': 'GeneXpert Dx',
    'assay_precision': 'Semi-Quantitative',
    'sdmc_pipeline_version':'Complex',
}

data_path1 = "/trials/covpn/p3008t/s001/qdata/LabData/MTB-Rif_pass-through/CoVPN3008TB_79004_2025.12.22_12.15.02.csv"
df1 = process_complex_submission(data_path1)

# Merge on guspecs
guspecs = pd.read_csv('/trials/covpn/p3008t/s001/qdata/LabData/MTB-Rif_pass-through/3008_Global specs.csv')
guspecs = guspecs.rename(columns={'Global Unique Id':'guspec'})
guspecs = guspecs.dropna(how='all')
guspecs.ptid = guspecs.ptid.astype(int).astype(str)

df1['ptid'] = df1['Sample ID'].str[:9]
df1.ptid = df1.ptid.astype(str)

tmp = df1.merge(guspecs[['guspec','ptid','visitno']], on='ptid', how='left')
assert len(tmp) == len(df1)
df1 = df1.merge(guspecs[['guspec','ptid','visitno']], on='ptid', how='left')

df1 = df1.rename(columns={
    'Assay Version':'assay_version',
    'Assay Type':'assay_categorization_lab',
})
df1 = df1.drop(columns=['Assay','Test Type'])

outputs1 = sdmc.standard_processing(
    input_data=df1,
    input_data_path=data_path1,
    guspec_col="guspec",
    network="CoVPN",
    metadata_dict=md,
    ldms=ldms
)

outputs1[['result_mtb','result_rif']].drop_duplicates()

data_path2 = '/trials/covpn/p3008t/s001/qdata/LabData/MTB-Rif_pass-through/CoVPN3008_TB Substudy_79001_2025.12.22_12.02.01.csv'
df2 = process_complex_submission(data_path2)

# Merge on guspecs
df2['ptid'] = df2['Sample ID'].str[:9]
df2.ptid = df2.ptid.astype(str)

tmp = df2.merge(guspecs[['guspec','ptid','visitno']], on='ptid', how='left')
assert len(tmp) == len(df2)
df2 = df2.merge(guspecs[['guspec','ptid','visitno']], on='ptid', how='left')

df2 = df2.rename(columns={
    'Assay Version':'assay_version',
    'Assay Type':'assay_categorization_lab',
})
df2 = df2.drop(columns=['Assay','Test Type'])

outputs2 = sdmc.standard_processing(
    input_data=df2,
    input_data_path=data_path2,
    guspec_col="guspec",
    network="HVTN",
    metadata_dict=md,
    ldms=ldms
)

# outputs2[['result_mtb','result_rif']].drop_duplicates()

# this one has two rows in ldms
# ldms.loc[ldms.guspec=="0410-0SYFJA00-001"]

outputs = pd.concat([outputs1, outputs2])

assert len(outputs.loc[(outputs.guspec.notna()) & (outputs.visitno_x!=outputs.visitno_y.astype(float)),['visitno_x','visitno_y']]) == 0

assert len(outputs.loc[(outputs.guspec.notna()) & (outputs.ptid_x.astype(str)!=outputs.ptid_y.astype(str)),['ptid_x','ptid_y']]) == 0

outputs = outputs.drop(columns=['ptid_x','visitno_x'])
outputs = outputs.rename(columns={'ptid_y':'ptid', 'visitno_y':'visitno'})

# ADD IN ADDITIONAL RENAMING / COLUMN DROPPING HERE
outputs = outputs.rename(columns={'s/w_version':'lab_software_version'})
outputs.loc[outputs.guspec.isna(),'specrole'] = "Control"

reorder = [
    'network',
    'protocol',
    'specrole',
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
    'assay_details',
    'buffer_diluent',
    'buffer_dilution',
    'lab_software',
    'lab_software_version',
    'instrument',
    'assay_precision',
    'sample_id',
    'assay_version',
    'assay_categorization_lab',
    'test_disclaimer',
    'result_mtb',
    'result_rif',
    'analyte_name',
    'ct',
    'endpt',
    'analyte_result',
    'probe_check_result',
    'cartridge_index',
    'module_name',
    'module_s/n',
    'cartridge_s/n',
    'instrument_s/n',
    'reagent_lot_id',
    'expiration_date',
    'start_time',
    'end_time',
    'error_status',
    'status',
    'user',
    'sdmc_pipeline_version',
    'sdmc_processing_datetime',
    'sdmc_data_receipt_datetime',
    'input_file_name',
]
assert set(reorder).symmetric_difference(outputs.columns) == set()

outputs = outputs[reorder]

today = datetime.date.today().isoformat()
savedir = "/networks/vtn/lab/SDMC_labscience/studies/CoVPN/CoVPN3008/assays/MTB-Rif/misc_files/data_processing/working_outputs/"
fname=f"MTB_Rif_complex_output_processed_{today}.txt"

outputs.to_csv(
    savedir + fname, sep='\t', index=False
)

# 4-3-2026 i'm updaing this code to fix the specrole column. this check shows that only specrole and processing date changed:
compare = pd.read_csv(
    '/networks/vtn/lab/SDMC_labscience/studies/CoVPN/CoVPN3008/assays/MTB-Rif/misc_files/data_processing/archive/DRAFT_MTB_Rif_complex_output_processed_2026-03-31.txt',
    sep="\t"
)

outputs.ptid = outputs.ptid.astype(float)
outputs.visitno = outputs.visitno.astype(float)
outputs.lab_software_version = outputs.lab_software_version.astype(float)

for col in ['assay_version','cartridge_index', 'module_s/n', 'cartridge_s/n', 'instrument_s/n', 'reagent_lot_id']:
    outputs[col] = outputs[col].astype(int)
    
aa = outputs.reset_index(drop=True).compare(compare)

print(aa)