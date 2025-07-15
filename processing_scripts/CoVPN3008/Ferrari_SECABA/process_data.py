## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 07/15/2025
# Purpose:  - Process SECABA data from Ferrari lab
## ---------------------------------------------------------------------------##
import pandas as pd
import numpy as np
import datetime, os
import yaml
import sdmc_tools.process as sdmc
import sdmc_tools.constants as constants
import sdmc_tools.access_ldms as access_ldms

def main():
    # read in data
    secaba_analysis_path = "/trials/covpn/p3008/s001/qdata/LabData/SECABA_pass-through/Ferrari_CoVPN3008_SECABA_Analysis_2024-10-11.csv"
    secaba_metadata_path = "/trials/covpn/p3008/s001/qdata/LabData/SECABA_pass-through/Ferrari_CoVPN3008_SECABA_Metadata_2024-09-23.csv"

    data = pd.read_csv(secaba_analysis_path)
    metadata = pd.read_csv(secaba_metadata_path)

    data.columns = [i.lower().replace(" ", "_").replace("%", "pct_") for i in data.columns]

    rename = {'sample_id':'guspec', 'ptid':'submitted_ptid', 'visit': 'submitted_visit'}
    data = data.rename(columns=rename)
    metadata = metadata.rename(columns=rename)

    to_drop = metadata.loc[(metadata.guspec=="0392-01JWBJ00-002") & (metadata.fcs_file_name.str.contains("347800695"))]
    metadata = metadata.loc[~metadata.index.isin(to_drop.index)]

    # merge metadata onto data
    data = data.merge(metadata, on=['guspec', 'submitted_ptid', 'submitted_visit', 'dilution'], how='outer')

    ## hand-entered metadata
    metadata_dict = {
        'network': 'CoVPN',
        'upload_lab_id': 'GF',
        'assay_lab_name': 'Ferrari Lab',
        'instrument': 'Fortessa Flow Cytometer',
        'assay_type': 'SECABA',
        'specrole': 'Sample',
        'result_units': 'Percent',
    }

    # pull in ldms
    ldms = access_ldms.pull_one_protocol('covpn', 3008)

    ldms = ldms.loc[ldms.guspec.isin(data.guspec)]
    ldms.lstudy = ldms.lstudy.astype(int)
    ldms.txtpid = ldms.txtpid.astype(str).astype(int)

    ldms = ldms.drop_duplicates()

    outputs = sdmc.standard_processing(
        input_data=data,
        input_data_path=secaba_analysis_path,
        guspec_col='guspec',
        network='covpn',
        metadata_dict=metadata_dict,
        ldms=ldms,
        additional_input_paths={'input_metadata_file_name': secaba_metadata_path}

    )
    outputs = outputs.drop(columns=["submitted_visit", "submitted_ptid"])

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
        'instrument',
        'dilution',
        'mock_pct_igg+',
        'transfected_pct_igg+',
        'background_subtracted_pct_igg+',
        'result_units',
        'isotype',
        'assay_date',
        'operator',
        'analysis_file_name',
        'fcs_file_name',
        'xml_file_name',
        'sdmc_processing_datetime',
        'sdmc_data_receipt_datetime',
        'input_file_name',
        'input_metadata_file_name',
    ]
    assert(len(set(reorder).symmetric_difference(outputs.columns))==0)
    outputs = outputs[reorder]

    # make sure matches previous format
    outputs.protocol = outputs.protocol.astype(int)
    outputs.ptid = outputs.ptid.astype(int)
    outputs.visitno = outputs.visitno.astype(float)

    # save to .txt
    savedir = "/networks/vtn/lab/SDMC_labscience/studies/CoVPN/CoVPN3008/assays/SECABA/misc_files/data_processing/"
    today = datetime.datetime.today().strftime('%Y-%m-%d')
    fname = f"CoVPN3008_Ferrari_SECABA_processed_{today}.txt"
    outputs.to_csv(savedir + fname, index=False, sep="\t")

    old = pd.read_csv(savedir + '/archive_2025_07_15/CoVPN3008_Ferrari_SECABA_processed_2024-10-11.txt', sep="\t")
    diff = outputs.compare(old)
    diff['visitno'].dropna()

    ## pivot summaries
    detail_summary = pd.pivot_table(
        outputs,
        index=['ptid','visitno'],
        columns=['isotype','dilution'],
        aggfunc='count',
        fill_value=0,
    )[['background_subtracted_pct_igg+']].droplevel(level=0, axis=1)
    detail_summary.to_excel(savedir + "CoVPN3008_SECABA_summary_ptid_visitno_dilution_isotype.xlsx")

    simple_summary = pd.pivot_table(
        outputs,
        index='ptid',
        columns='visitno',
        aggfunc='count',
        fill_value=0
    )[['background_subtracted_pct_igg+']]
    simple_summary.to_excel(savedir + "CoVPN3008_SECABA_summary_ptid_visitno.xlsx")

    # check background subtraction
    check = outputs['transfected_pct_igg+'] - outputs['mock_pct_igg+']
    check[check < 0] = 0

    np.abs(check - outputs['background_subtracted_pct_igg+']).max()

    # check = outputs['Transfected %IgG+'] - data1011['Mock %IgG+']
    # check[check < 0] = 0

    # np.abs(check - data1011['Background subtracted %IgG+']).max()

    ## additional checks
    manifest_path = "/networks/vtn/lab/SDMC_labscience/studies/CoVPN/CoVPN3008/assays/SECABA/misc_files/CoVPN_shipping_manifest.txt"
    manifest = pd.read_csv(manifest_path, sep="\t")

    manifest_path2 = "/networks/vtn/lab/SDMC_labscience/studies/CoVPN/CoVPN3008/assays/SECABA/misc_files/CoVPN_shipping_manifest_pt2.txt"
    manifest2 = pd.read_csv(manifest_path2, sep="\t")

    manifest_ids = set(manifest2.loc[manifest2.PROTOCOL==3008.].GLOBAL_ID).union(manifest.GLOBAL_ID)
    check = manifest_ids.symmetric_difference(outputs.guspec)

    assert(len(check)==0)

    # len(set(outputs.guspec).difference(manifest.GLOBAL_ID))
    # outputs.guspec.nunique()
    # outputs.ptid.nunique()
    # outputs.groupby(['ptid']).visitno.nunique().value_counts()
    # manifest.GLOBAL_ID.nunique()

    # for_nick = outputs[['guspec','ptid','visitno']].drop_duplicates()
    # for_nick["in_manifest"] = for_nick.guspec.isin(manifest.GLOBAL_ID)
    # for_nick.in_manifest.sum()

    # savepath = "/networks/vtn/lab/SDMC_labscience/operations/documents/templates/assay/template_testing/CoVPN3008_Ferrari_samples.txt"
    # for_nick.to_csv(savepath, sep="\t", index=False)

    # # visitno per samples in manifest
    # outputs.loc[outputs.guspec.isin(manifest.GLOBAL_ID)].visitno.value_counts()

    # # visitno per samples not in manifest
    # outputs.loc[~outputs.guspec.isin(manifest.GLOBAL_ID)].visitno.value_counts()

    # outputs.loc[outputs.guspec.isin(list(set(outputs.guspec).difference(manifest.GLOBAL_ID)))]

    # data1011 = pd.read_csv("/trials/covpn/p3008/s001/qdata/LabData/SECABA_pass-through/Ferrari_CoVPN3008_SECABA_Analysis_2024-10-11.csv")
    # metadata0923 = pd.read_csv("/trials/covpn/p3008/s001/qdata/LabData/SECABA_pass-through/Ferrari_CoVPN3008_SECABA_Metadata_2024-09-23.csv")

    # data0923 = pd.read_csv("/trials/covpn/p3008/s001/qdata/LabData/SECABA_pass-through/archive/Ferrari_CoVPN3008_SECABA_Analysis_2024-09-23.csv")

    # data0102 = pd.read_csv("/trials/covpn/p3008/s001/qdata/LabData/SECABA_pass-through/quarantine/Ferrari_CoVPN3008_SECABA_Analysis_2024-01-02.csv")
    # data0612 = pd.read_csv("/trials/covpn/p3008/s001/qdata/LabData/SECABA_pass-through/quarantine/Ferrari_CoVPN3008_SECABA_Analysis_2024-06-12.csv")

    # data0102.compare(data0612)
    # data0612.compare(data0923)
    # data0923.compare(data1011)


    # check = data1011['Transfected %IgG+'] - data1011['Mock %IgG+']
    # check[check < 0] = 0

    # np.abs(check - data1011['Background subtracted %IgG+']).max()

if __name__ == '__main__':
    main()
