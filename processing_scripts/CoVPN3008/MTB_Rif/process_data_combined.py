## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 04/03/2026
# Purpose: Process MTB Rif simple machine outputs/converted pdfs 
##         and merge on complex machine outputs
## ---------------------------------------------------------------------------##
import pandas as pd
import numpy as np
import sdmc_tools.access_ldms as access_ldms
import sdmc_tools.process as sdmc
import os
import datetime

from io import StringIO

import pandas as pd
import numpy as np
import sdmc_tools.access_ldms as access_ldms
import sdmc_tools.process as sdmc
import os
import datetime

from io import StringIO

def main():
    ## read in data
    pdf_dir = '/trials/covpn/p3008t/s001/qdata/LabData/MTB-Rif_pass-through/'
    pdf_data = [
        '784801133V302.xlsx',
        '792801574V302.xlsx',
        '807803673V302.xlsx',
        '816802568V302.xlsx',
        '836803289V302.xlsx',
        '836803582V302.xlsx',
    ]
    
    pdfs = pd.concat([parse_pdf(pdf_dir + p) for p in pdf_data])
    
    ## read in simple outputs ---------------------------------------------------------------------------- ##
    datadir = '/trials/covpn/p3008t/s001/qdata/LabData/MTB-Rif_pass-through/'
    
    input_data_path1 = '/trials/covpn/p3008t/s001/qdata/LabData/MTB-Rif_pass-through/66_Tests_Exported_2026.01.20_13.21.27.csv'
    df1 = pd.read_csv(input_data_path1)
    
    input_data_path2 = '/trials/covpn/p3008t/s001/qdata/LabData/MTB-Rif_pass-through/68_Tests_Exported_2026.01.20_13.29.00.csv'
    df2 = pd.read_csv(input_data_path2)
    
    # save which rows map to which file
    input_file_map = {i:'66_Tests_Exported_2026.01.20_13.21.27.csv' for i in df1['Patient ID']}
    input_file_map = input_file_map | {i:'68_Tests_Exported_2026.01.20_13.29.00.csv' for i in df2['Patient ID']}
    
    sdmc_data_receipt_time_map1 = {i:datetime.datetime.fromtimestamp(os.path.getmtime(input_data_path1)).isoformat() for i in df1['Patient ID']}
    sdmc_data_receipt_time_map2 = {i:datetime.datetime.fromtimestamp(os.path.getmtime(input_data_path2)).isoformat() for i in df2['Patient ID']}
    sdmc_data_receipt_time_map = sdmc_data_receipt_time_map1 | sdmc_data_receipt_time_map2
    
    df1 = df1.dropna(how='all')
    df2 = df2.dropna(how='all')
    
    df = pd.concat([df1, df2])
    df = pd.concat([
        df.drop(columns='Result'),
        df.Result.str.split("; ", expand=True).rename(columns={0:'result_mtb_interpretation',1:'result_rif_interpretation'})
    ], axis=1)
    df.columns = [i.lower().replace(" ","_") for i in df.columns]
    df = df.rename(columns={
        'reagent_lot':'reagent_lot_id',
        'sample_id':'lab_internal_sample_id_date',
        'patient_id':'sample_id'
    })
    df = df.drop(columns=['assay_name','test_type'])
    
    df.reagent_lot_id = df.reagent_lot_id.astype(float)
    pdfs.reagent_lot_id = pdfs.reagent_lot_id.astype(float)
    
    set(df.columns).intersection(pdfs.columns)
    
    tmp = df.merge(pdfs, on=['sample_id','reagent_lot_id'], how='left')
    
    # THE MERGE SHOULD ADD ROWS FOR THE PTIDS IN THE PDFS
    assert len(df) - pdfs.sample_id.nunique() + len(pdfs) == len(tmp)
    df = df.merge(pdfs, on=['sample_id','reagent_lot_id'], how='left')
        
    # Merge on guspecs
    guspecs = pd.read_csv(datadir + '3008_Global specs.csv')
    guspecs = guspecs.rename(columns={'Global Unique Id':'guspec'})
    guspecs = guspecs.dropna(how='all')
    
    df['ptid'] = df.sample_id.str[:-4]
    df.ptid = df.ptid.astype(int)
    guspecs.ptid = guspecs.ptid.astype(int)
    
    tmp = df.merge(guspecs[['guspec','ptid','visitno']], on='ptid', how='left')
    assert len(tmp) == len(df)
    df = df.merge(guspecs[['guspec','ptid','visitno']], on='ptid', how='left')
    
    ## standard processing ------------------------------------------------------------------------------- ##
    ldms = access_ldms.pull_one_protocol('covpn', 3008)
    
    md = {
        'network': 'CoVPN',
        'protocol': 3008,
        'specrole':'Sample',
        'upload_lab_id': 'CH',
        'assay_lab_name': 'CHIL',
        'assay_type': 'RT-PCR', 
        'assay_subtype': 'Xpert MTB/RIF Ultra',
        'assay_details':'Assay modified for tongue swab specimens by diluting SR Buffer',
        'buffer_diluent':'Phosphate Buffered Saline (PBS)',
        'buffer_dilution':'66%',
        'instrument': 'P049 Xpert MTB/RIF Ultra', 
        'lab_software': 'GeneXpert Dx',
        'assay_precision': 'Semi-Quantitative',
        'sdmc_pipeline_version':'Simple',
    }
    
    # @sara the real guspecs should be merged on here!
    # @sara in addition the real filepaths to the input data need to be merged onto each row - im sorry i didnt get to this
    outputs = sdmc.standard_processing(
        input_data=df,
        input_data_path='FILL THIS IN',
        guspec_col='guspec',
        network="CoVPN",
        metadata_dict=md,
        ldms=ldms
    )
    
    set(outputs.sample_id).difference(input_file_map.keys())
    
    outputs['input_file_name'] = outputs.sample_id.map(input_file_map)
    outputs['sdmc_data_receipt_datetime'] = outputs.sample_id.map(sdmc_data_receipt_time_map)
    
    outputs.ptid_x = outputs.ptid_x.astype(int)
    outputs.ptid_y = outputs.ptid_y.astype(int)
    assert (outputs.ptid_x!=outputs.ptid_y).sum() == 0
    
    outputs.visitno_x = outputs.visitno_x.astype(float)
    outputs.visitno_y = outputs.visitno_y.astype(float)
    assert (outputs.visitno_x!=outputs.visitno_y).sum() == 1 #0332-02PJVB00-001 has a weird value in the input file
    
    ldms_discrep = outputs.loc[outputs.visitno_x!=outputs.visitno_y,['guspec','visitno_x','visitno_y']]
    ldms_discrep = ldms_discrep.rename(columns={'visitno_x':'visitno_ref_sheet', 'visitno_y':'visitno_ldms'})
    
    outputs = outputs.drop(columns=['ptid_x','visitno_x'])
    outputs = outputs.rename(columns={'ptid_y':'ptid', 'visitno_y':'visitno'})
    
    ## merge on complex outputs ---------------------------------------------------------------------------- ##
    complex_outputs = '/networks/vtn/lab/SDMC_labscience/studies/CoVPN/CoVPN3008/assays/MTB-Rif/misc_files/data_processing/working_outputs/MTB_Rif_complex_output_processed_2026-04-03.txt'
    complex_outputs = pd.read_csv(complex_outputs, sep='\t')
    complex_outputs = complex_outputs.drop(columns=['assay_categorization_lab'])
    
    complex_outputs = complex_outputs.rename(columns={
        'analyte_result':'result_qualitative',
        'result_mtb':'result_mtb_interpretation',
        'result_rif':'result_rif_interpretation',
        'cartridge_s/n':'cartridge_serial',
        'module_s/n':'module_serial',
        'instrument_s/n':'instrument_serial',
        'probe_check_result':'result_probe_control',
        'expiration_date':'reagent_lot_expiration',
        'endpt':'result_fi',
    })
    
    # these contain almost the same columns -- no cartridge_index in the simple machine outputs
    # set(reorder).difference(complex_outputs.columns)
    # set(complex_outputs.columns).difference(reorder)
    # set(complex_outputs.columns).difference(outputs.columns)
    
    outputs = pd.concat([complex_outputs, outputs])
    outputs = outputs.rename(columns={
        'analyte_name':'target',
        'reagent_lot_id':'reagent_lot',
        'start_date_&_time':'start_datetime',
        'start_time':'test_datetime',
        'sample_id':'lab_internal_sample_id',
    })
    outputs = outputs.drop(columns='end_time')
    
    outputs['qc_flag_sdmc'] = False
    outputs.loc[outputs.guspec=="0410-0SYFJA00-001", 'qc_flag_sdmc'] = True
    outputs.loc[outputs.guspec=="0410-0SYFJA00-001", 'qc_flag_sdmc_detail'] = "This guspec has two different drawdt values in LDMS. Both versions have been included here, so this sample has 12 records instead of 6. We don't know which is correct and are looking into this."
    
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
        'assay_precision',
        'instrument',
        'instrument_serial',
        'lab_software',
        'lab_software_version',
        'target',
        'assay_version',
        'cartridge_index',
        'cartridge_serial',
        'ct',
        'error_status',
        'lab_internal_sample_id_date',
        'module_name',
        'module_serial',
        'reagent_lot',
        'reagent_lot_expiration',
        'result_qualitative',
        'result_mtb_interpretation',
        'result_rif_interpretation',
        'result_probe_control',
        'result_fi',
        'lab_internal_sample_id',
        'start_datetime',
        'test_datetime',
        'status',
        'test_disclaimer',
        'user',
        'sdmc_pipeline_version',
        'qc_flag_sdmc',
        'qc_flag_sdmc_detail',
        'sdmc_processing_datetime',
        'sdmc_data_receipt_datetime',
        'input_file_name',
        'input_file_converted_pdf',
    ]
    
    set(reorder).symmetric_difference(outputs.columns)
    
    assert set(reorder).symmetric_difference(outputs.columns) == set()
    outputs = outputs[reorder]
    outputs.loc[outputs.guspec.isna(), 'specrole'] = 'Control'
    outputs.result_qualitative = outputs.result_qualitative.replace(["nan","NA"], np.nan)
    for col in ['instrument_serial','cartridge_serial','ct','module_serial','result_fi','ptid','visitno','lab_software_version']:
        outputs[col] = outputs[col].astype(float)
    
    
    usecols = pd.read_excel("/networks/vtn/lab/SDMC_labscience/assays/Xpert MTB_RIF Ultra/SDMC_materials/DRAFT_qdata_format_MTB-Rif.xlsx")
    set(usecols.variable_name).difference(outputs.columns)
    set(reorder).difference(usecols.variable_name)
    
    # DIF CHECK TO SHOW THAT NOTHING CHANGED FROM THE MARCH 31 VERSION OTHER THAN LAB_INTERNAL_SAMPLE_ID
    outputs_03_31 = pd.read_csv(
        '/networks/vtn/lab/SDMC_labscience/studies/CoVPN/CoVPN3008/assays/MTB-Rif/misc_files/data_processing/archive/CoVPN3008_MTB-Rif_combined_output_processed_2026-03-31.txt',
        sep="\t"
    )
    
    assert set(outputs.columns).symmetric_difference(outputs_03_31.columns) == {'lab_internal_sample_id', 'sample_id'}
    
    outputs_03_31 = outputs_03_31.rename(columns={'sample_id':'lab_internal_sample_id'})
    outputs.ptid = outputs.ptid.astype(float)
    outputs.visitno = outputs.visitno.astype(float)
    outputs.lab_software_version = outputs.lab_software_version.astype(float)
            
    compare = outputs.reset_index(drop=True).compare(outputs_03_31.reset_index(drop=True))
    
    print(compare)
    
    ## save ----------------------------------------------------------------------------------------------- ##
    today = datetime.date.today().isoformat()
    savedir = "/networks/vtn/lab/SDMC_labscience/studies/CoVPN/CoVPN3008/assays/MTB-Rif/misc_files/data_processing/"
    fname=f"MTB_Rif_combined_output_processed_{today}.txt"
    
    outputs.to_csv(
        savedir + fname, sep='\t', index=False
    )
    
    # pivot_summary = pd.pivot_table(
    #     outputs,
    #     index='ptid',
    #     columns='visitno',
    #     values='result_mtb_interpretation',
    #     aggfunc='count',
    #     dropna=False
    # ).fillna(0)
    
    # ldms_discrep = ldms_discrep.merge(outputs.loc[outputs.guspec=="0332-02PJVB00-001",['ptid','drawdt','guspec']], on='guspec')
    
    # with pd.ExcelWriter(savedir + f'CoVPN3008_MTB_Rif_pivot_summary_{today}.xlsx', engine="openpyxl") as writer:
    #     pivot_summary.to_excel(writer, sheet_name="pivot_summary", index=True)
    #     ldms_discrep.to_excel(writer, sheet_name="ldms_discrep", index=True)
    


    # manifest check 
    manifest = pd.read_csv(
        '/networks/vtn/lab/SDMC_labscience/studies/CoVPN/CoVPN3008/assays/MTB-Rif/FH_shipping_manifest.csv',
    )
    
    set(outputs.guspec).difference(manifest['Global Spec ID'])
    
    len(set(manifest['Global Spec ID']).difference(outputs.guspec))
    
    outputs.guspec.nunique()
    
    manifest['guspec_core'] = manifest['Global Spec ID'].str.rpartition("-")[0]
    
    outputs['guspec_core'] = outputs.guspec.str.rpartition("-")[0]
    
    set(outputs.guspec_core).symmetric_difference(manifest.guspec_core)

# helper fn for rows of data that got combined into one excel cell separated by a "\n"
def split_row(row):
    # Split every column on newline
    split_cols = row.astype(str).str.split("\n")
    
    # Determine how many rows this will expand into
    max_len = split_cols.str.len().max()
    
    # Pad shorter splits with None so lengths match
    split_cols = split_cols.apply(lambda x: x + [None]*(max_len - len(x)))
    
    # Convert to DataFrame
    return pd.DataFrame(dict(zip(row.index, split_cols)))

## parse pdf data ---------------------------------------------------------------------------- ##
def parse_pdf(pdf_path):
    data = pd.read_excel(
        pdf_path,
        header=None
    )
    
    metadata = {}
    blocks = []
    
    metadata['sample_id'] = data.iloc[0].dropna().tolist()[1].split("\n")[0]
    
    tmp = data.iloc[2:4].dropna(how='all', axis=1)
    metadata_new = dict(zip(tmp.iloc[0], tmp.iloc[1]))
    metadata = metadata | metadata_new
        
    data.loc[data[0]=="Reagent Lot*: Error Status:", 0] = "Reagent Lot*:\nError Status:"
    idx_start = data[(
        (data.apply(lambda row: row.astype(str).str.replace(" ","").str.lower().str.contains("user:").any(), axis=1)) &
        (data.apply(lambda row: row.astype(str).str.replace(" ","").str.lower().str.contains("status").any(), axis=1))
    )].index[0]
    idx_end = data[data.apply(lambda row: row.astype(str).str.replace(" ","").str.lower().str.contains("reagentlot").any(), axis=1)].index[0]
    tmp = data.iloc[idx_start:idx_end+1].dropna(how='all', axis=1)
    tmp = pd.concat([split_row(row) for _, row in tmp.iterrows()],
                  ignore_index=True)
    tmp.columns = ['a','b','a','b']
    tmp = pd.concat([
        tmp.iloc[:,[0,1]],
        tmp.iloc[:,[2,3]]
    ]).dropna(how='all')
    metadata_new = dict(zip(tmp.iloc[:,0], tmp.iloc[:,1]))
    metadata = metadata | metadata_new
           
    idx_start = data[data.apply(lambda row: row.astype(str).str.replace(" ","").str.lower().str.contains("analytect").any(), axis=1)].index[0]
    idx_end = data[(
        (data.apply(lambda row: row.astype(str).str.replace(" ","").str.lower().str.contains("user:").any(), axis=1)) &
        (data.apply(lambda row: row.astype(str).str.replace(" ","").str.lower().str.contains("status").any(), axis=1))
    )].index[0]
    df = data.iloc[idx_start+1:idx_end].dropna(how='all', axis=1)
    df.columns = ['analyte','ct','endpoint','analyte_result','result_probe_control']
    df = pd.concat([split_row(row) for _, row in df.iterrows()],
                  ignore_index=True)
    
    df = df.merge(pd.DataFrame([metadata]), how='cross')
    disclaimer = data[data.apply(lambda row: row.astype(str).str.replace(" ","").str.lower().str.contains("diagnostic").any(), axis=1)].values[-1][0]
    df['test_disclaimer'] = disclaimer
    df.columns = [i.replace("*",'') for i in df.columns]
    column_rename = {
        'analyte': 'analyte_name',
        # 'ct': 'ct',
        'endpoint': 'result_fi',
        'analyte_result': 'result_qualitative',
        # 'result_probe_control': 'result_probe_control',
        # 'sample_id': 'sample_id',
        # 'Assay': 'Assay',
        'Assay Name':'assay_name',
        'Assay Version': 'assay_version',
        'Assay Type': 'assay_categorization_lab',
        'User: Status:': 'user_status',
        'User:':'user',
        'Status:':'status',
        'Error Status:':'error_status',
        'Reagent Lot ID: Notes:':'reagent_lot_id',
        'Reagent Lot:':'reagent_lot_id',
        'Expiration Date:': 'reagent_lot_expiration',
        'S/W Version:': 'lab_software_version',
        'Cartridge S/N:': 'cartridge_serial',
        'Start Time:': 'start_time',
        'End Time:': 'end_time',
        'Instrument S/N:': 'instrument_serial',
        'Module S/N:': 'module_serial',
        'Module Name:':'module_name',
    }
    df = df.rename(columns=column_rename)
    df.columns = [i.lower().replace(" ","_") for i in df.columns]
    dropcols = [
        'assay', 'assay_name','test_type', 'assay_categorization_lab'
    ]
    dropcols = list(set(dropcols).intersection(df.columns))
    df = df.drop(columns=dropcols)
    if 'user_status' in df.columns:
        df = pd.concat([
            df.drop(columns='user_status'),
            df.user_status.str.split(" ", expand=True).rename(columns={0:'user',1:'status'})
        ], axis=1)
    df['input_file_converted_pdf'] = pdf_path.split("/")[-1]
    return df

if __name__=="__main__":
    main()