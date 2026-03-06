import pandas as pd
import numpy as np
import sdmc_tools.access_ldms as access_ldms
import sdmc_tools.process as sdmc
import os
import datetime

from io import StringIO


def main():
    ## parse pdf data ---------------------------------------------------------------------------- ##
    def parse_pdf(pdf_path):
        data = pd.read_excel(
            pdf_path,
            header=None
        )
        
        metadata = {}
        blocks = []
        # blocks += [data.iloc[0].dropna().tolist()]
        
        metadata['sample_id'] = data.iloc[0].dropna().tolist()[1].split("\n")[0]
        
        # blocks += [data.iloc[2:4].dropna(how='all', axis=1)]
        tmp = data.iloc[2:4].dropna(how='all', axis=1)
        metadata_new = dict(zip(tmp.iloc[0], tmp.iloc[1]))
        metadata = metadata | metadata_new
        
        # blocks += [data.iloc[5].dropna().tolist()]
        # metadata['result_mtb_interpretation'] = data.iloc[5].dropna().tolist()[1]
        
        idx_start = data[data.apply(lambda row: row.astype(str).str.replace(" ","").str.lower().str.contains("user:status").any(), axis=1)].index[0]
        idx_end = data[data.apply(lambda row: row.astype(str).str.replace(" ","").str.lower().str.contains("reagentlotid").any(), axis=1)].index[0]
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
        idx_end = data[data.apply(lambda row: row.astype(str).str.replace(" ","").str.lower().str.contains("user:status:").any(), axis=1)].index[0]
        df = data.iloc[idx_start+1:idx_end].dropna(how='all', axis=1)
        df.columns = ['analyte','ct','endpoint','analyte_result','result_probe_control']
        df = pd.concat([split_row(row) for _, row in df.iterrows()],
                      ignore_index=True)
        df = df.merge(pd.DataFrame([metadata]), how='cross')
        disclaimer = data[data.apply(lambda row: row.astype(str).str.replace(" ","").str.lower().str.contains("diagnostic").any(), axis=1)].values[-1][0]
        df['test_disclaimer'] = disclaimer
        column_rename = {
            'analyte': 'analyte_name',
            # 'ct': 'ct',
            'endpoint': 'endpt',
            'analyte_result': 'result_qualitative',
            # 'result_probe_control': 'result_probe_control',
            # 'sample_id': 'sample_id',
            # 'Assay': 'Assay',
            'Assay Version': 'assay_version',
            'Assay Type': 'assay_categorization_lab',
            'User: Status:': 'user_status',
            'Error Status:':'error_status',
            'Reagent Lot ID*: Notes:':'reagent_lot_id',
            'Expiration Date*:': 'reagent_lot_expiration',
            'S/W Version:': 'lab_software_version',
            'Cartridge S/N*:': 'cartridge_serial',
            'Start Time:': 'start_time',
            'End Time:': 'end_time',
            'Instrument S/N:': 'instrument_serial',
            'Module S/N:': 'module_serial',
            'Module Name:':'module_name',
        }
        df = df.rename(columns=column_rename)
        df.columns = [i.lower().replace(" ","_") for i in df.columns]
        dropcols = [
            'assay', 'assay_name','test_type'
        ]
        dropcols = list(set(dropcols).intersection(df.columns))
        df = df.drop(columns=dropcols)
        df = pd.concat([
            df.drop(columns='user_status'),
            df.user_status.str.split(" ", expand=True).rename(columns={0:'user',1:'status'})
        ], axis=1)
        return df

    ## @sara this is where the real pdf filepaths should go!
    pdf_dir = '/networks/vtn/lab/SDMC_labscience/assays/MTB-Rif/example_data/'
    pdf_data = [
        'pdf_conversion_123456789V302_2025.12.22_10.26.49_details.xlsx',
        'pdf_conversion_297800187V302_2025.12.22_10.26.49_details.xlsx',
        'pdf_conversion_339800837V302_2025.12.22_12.01.19_details.xlsx',
    ]

    pdfs = pd.concat([parse_pdf(pdf_dir + p) for p in pdf_data])

    ## @sara below i mapped the pdf ids to ones that are actually in the data to make the fake data look real.
    ## you should delete this and then check that the merge worked well -- ill @ you below
    # PRETENDING THESE ARE REAL IDS THAT WE WANTED! THIS IS CREATING FAKE DATA
    # THE REAL DATA SHOULD ALREADY HAVE MATCHING IDS
    pdfs = pdfs.loc[pdfs.sample_id!='339800837V302']
    pdfs.sample_id = pdfs.sample_id.map({
        '123456789V302':'792801574V302', 
        '297800187V302':'784801133V302',
    })

    ## read in simple outputs ---------------------------------------------------------------------------- ##
    datadir = '/trials/covpn/p3008t/s001/qdata/LabData/MTB-Rif_pass-through/'

    input_data_path1 = '/trials/covpn/p3008t/s001/qdata/LabData/MTB-Rif_pass-through/66_Tests_Exported_2026.01.20_13.21.27.csv'
    df1 = pd.read_csv(input_data_path1)

    input_data_path2 = '/trials/covpn/p3008t/s001/qdata/LabData/MTB-Rif_pass-through/68_Tests_Exported_2026.01.20_13.29.00.csv'
    df2 = pd.read_csv(input_data_path2)

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
    df = df.drop(columns=['assay_name'])

    df.reagent_lot_id = df.reagent_lot_id.astype(float)
    pdfs.reagent_lot_id = pdfs.reagent_lot_id.astype(float)
    tmp = df.merge(pdfs, on=['sample_id','reagent_lot_id'], how='left')

    # @sara ive just added this simple check to ensure the merge worked correctly. might want to check this more closely
    # THE MERGE SHOULD ADD ROWS FOR THE PTIDS IN THE PDFS
    assert len(df) - pdfs.sample_id.nunique() + len(pdfs) == len(tmp)
    df = df.merge(pdfs, on=['sample_id','reagent_lot_id'], how='left')

    ## standard processing ------------------------------------------------------------------------------- ##
    ldms = access_ldms.pull_one_protocol('covpn', 3008)

    md = {
        'network': 'CoVPN',
        'protocol': 3008,
        'specrole':'Sample',
        'upload_lab_id': 'CH',
        'assay_lab_name': 'CHIL',
        'assay_type': 'RT-PCR', 
        'assay_subtype': 'MTB/Rif',
        'assay_kit': 'Xpert MTB/RIF Ultra',
        'instrument': 'P049 Xpert MTB/RIF Ultra', 
        'lab_software': 'GeneXpert Dx',
        'assay_precision': 'Semi-Quantitative',
        'pipeline_version':'Complex',
    }

    # @sara the real guspecs should be merged on here!
    # @sara in addition the real filepaths to the input data need to be merged onto each row - im sorry i didnt get to this
    df['FAKE_PLACEHOLDER_GUSPEC'] = '0666-004MZG00-001'
    outputs = sdmc.standard_processing(
        input_data=df,
        input_data_path='FILL THIS IN',
        guspec_col="FAKE_PLACEHOLDER_GUSPEC",
        network="CoVPN",
        metadata_dict=md,
        ldms=ldms
    )

    ## merge on complex outputs ---------------------------------------------------------------------------- ##

    complex_outputs = '/networks/vtn/lab/SDMC_labscience/studies/CoVPN/CoVPN3008/assays/MTB-Rif/misc_files/data_processing/DRAFT_MTB_Rif_complex_output_processed_2026-03-04.txt'
    complex_outputs = pd.read_csv(complex_outputs, sep='\t')
    complex_outputs['pipeline_version'] = 'Simple'
    complex_outputs['network'] = 'CoVPN'
    outputs = pd.concat([complex_outputs, outputs.drop(columns='fake_placeholder_guspec')])

    ## save ----------------------------------------------------------------------------------------------- ##
    today = datetime.date.today().isoformat()
    savedir = "/networks/vtn/lab/SDMC_labscience/studies/CoVPN/CoVPN3008/assays/MTB-Rif/misc_files/data_processing/"
    fname=f"DRAFT_FAKE_MTB_Rif_combined_output_processed_{today}.txt"

    outputs.to_csv(
        savedir + fname, sep='\t', index=False
    )


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


if __name__=="__main__":
    main()