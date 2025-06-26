import pandas as pd
import numpy as np
import datetime
import yaml
import os
import sdmc_tools.constants as constants
import smtplib
from sdmc_tools.access_ldms import pull_one_protocol as get_ldms
from email.message import EmailMessage

# ---------------------------------------------------------------------------- #
yamlpath = '/home/bhaddock/repos/sdmc-adhoc/constants.yaml'
with open(yamlpath, 'r') as file:
    yamldict = yaml.safe_load(file)

def main():
    # this is everything since BH joined (jan 2024)
    endpoints = find_endpoints("/home/bhaddock/repos/sdmc-adhoc/processing_scripts", l=[])
    yams = [i for i in endpoints if i[-4:]=="yaml"]
    output_paths = [get_output_path(i) for i in yams]

    # these are the outputs known to not contain a guspec column; cant check these against ldms
    no_guspec = [
        '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN302/assays/mAb_ICABA/misc_files/data_processing/DRAFT_HVTN302_Ferrari_mAb_ICABA_processed_2024-07-30.txt',
        '/networks/vtn/lab/SDMC_labscience/studies/VISC/G002/MRGPRX2/data_processing/CAVD_G002_MRGPRX2_Processed_2025-02-27.txt',
        '/networks/cavd/VDCs/Schief/Schief_856-G002/SkinReactions/data/Glycan_array_Scripps/processed_data/CAVD_G002_Glycan_Microarray_data_processed_2025-01-23.txt',
        '/networks/cavd/VDCs/Schief/Schief_856-G002/SkinReactions/data/Glycan_array_Scripps/processed_data/CAVD_G002_Glycan_Microarray_data_processed_2025-01-23.txt',
        '/networks/vtn/lab/SDMC_labscience/studies/VISC/ID52-Frahm/omics_misc/misc_files/data_processing/GCMS_Metabolomics/GCMS_Metabolomics_processed_2025-03-04.txt',
        '/networks/vtn/lab/SDMC_labscience/studies/VISC/ID52-Frahm/omics_misc/misc_files/data_processing/LCMS_Metabolomics/LCMS_Metabolomics_processed_2025-03-04.txt',
        '/networks/vtn/lab/SDMC_labscience/studies/VISC/ID52-Frahm/omics_misc/misc_files/data_processing/LCMS_Proteomics/LCMS_Proteomics_processed_2025-03-11.txt',
        '/networks/vtn/lab/SDMC_labscience/studies/VISC/ID52-Frahm/omics_misc/misc_files/data_processing/LCMS_Lipidomics/LCMS_Lipidomics_processed_2025-03-04.txt'
        '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN137/assays/ICABA/misc_files/data_processing_mAb/HVTN137_mAb_ICABA_processed_2025-06-25.txt'
    ]

    # these are all the pre-2024 files
    olds = [
        '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN Auxiliary Studies/111-ANS_Maenetje_HVTN_503/draft_FCM08_output_20210616.txt',
        '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN Auxiliary Studies/Mycoplasma/161-EXS_Mycoplasma_Mayo_HVTN_092_096/sdmc_output/161EXS_Mycoplasma_Mayo_HVTN_092_096_20190312.txt',
        '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN Auxiliary Studies/Mycoplasma/161-EXS_Mycoplasma_Mayo_HVTN_092_096/sdmc_output/161EXS_Mycoplasma_Mayo_HVTN_092_096_20190313.txt',
        '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN Auxiliary Studies/Mycoplasma/162-EXS_Mycoplasma_Iowa_HVTN_092_096/sdmc_output/162EXS_Mycoplasma_Iowa_HVTN_092_096_Part1_20190225.txt',
        '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN Auxiliary Studies/Mycoplasma/162-EXS_Mycoplasma_Iowa_HVTN_092_096/sdmc_output/162EXS_Mycoplasma_Iowa_HVTN_092_096_Part3_20190711.txt',
        '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN Auxiliary Studies/Mycoplasma/162-EXS_Mycoplasma_Iowa_HVTN_092_096/sdmc_output/162EXS_Mycoplasma_Iowa_HVTN_092_096_Part3_20191014.txt',
        '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN139/assays/NAb/Wistar/misc_files/HVTN139_Adenovirus_nAb_Wistar_2023-08-30.csv',
        '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN139/assays/NAb/Wistar/misc_files/archive/HVTN139_Adenovirus_nAb_Wistar_2023-07-17.csv',
        '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN302/assays/AE_assays/basophil_activation_testing/misc_files/HVTN302_basophil_activation_processed_2023-11-22.txt',
        '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN303/assays/AE_assays/Circulating Immune Complex Panel/misc_files/data_processing/DRAFT_HVTN303_ARUP_C1C-C1q_processed_2023-09-25.txt',
        '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN303/assays/MSD/Binding_antibody_via_MSD/misc_files/data_processing/DRAFT_HVTN303_MSD_binding_processed_2023-09-26.txt',
        '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN303/assays/MSD/Cytokines_via_MSD/misc_files/data_processing/HVTN303_MSD_cytokines_processed_2024-02-16.txt',
    ]

    # for each dataset, check against ldms
    results = {}
    for o in set(output_paths).difference(no_guspec).union(olds):
        try:
            results[o] = check_for_ldms_changes(o)
        except Exception as error:
            results[o] = 'encountered error'

    # collect any issues
    mismatches = [i for i in results if results[i]=='mismatch']
    missing_guspecs = [i for i in results if results[i]=='guspec/guspec_core col missing']
    errors = [i for i in results if results[i]=='encountered error']

    # write myself an email
    email_text = ''

    if len(mismatches) > 0:
        email_text += "Misalignment between the following datasets and LDMS:\n"
        for m in mismatches:
            email_text += "\t- " + m + "\n"

    if len(missing_guspecs) > 0:
        if len(email_text)>0:
            email_text += "\n"
        email_text += "The following datasets don't appear to have a guspec col. Please examine:\n"
        for m in missing_guspecs:
            email_text += "\t- " + m + "\n"

    if len(errors) > 0:
        if len(email_text)>0:
            email_text += "\n"
        email_text += "The code threw an error trying to examine the following datasets. Couldn't compare against LDMS; please examine:\n"
        for m in errors:
            email_text += "\t- " + m + "\n"


    email_subject = "LDMS Monitoring code requires attention"
    if len(mismatches) > 0:
        email_subject = "ALERT: dataset misalignment with LDMS"


    recipients = [
        'bhaddock@fredhutch.org',
        'lbunts@fredhutch.org',
        'sthiebau@fredhutch.org'
    ]

    # if there's anything to report on, send email
    if len(email_text) > 0:
        msg = EmailMessage()
        msg.set_content(email_text)

        msg['Subject'] = email_subject
        msg['From'] = 'bhaddock@fredhutch.org'
        msg['To'] = ", ".join(recipients)

        s = smtplib.SMTP('localhost')
        s.send_message(msg)
        s.quit()

# ---------------------------------------------------------------------------- #
def check_for_ldms_changes(path):
    if '.txt' in path:
        df = pd.read_csv(path, sep="\t")
    else:
        df = pd.read_csv(path)
    df.columns = [i.lower() for i in df.columns]
    if 'specrole' in df.columns:
        df = df.loc[df.specrole=="Sample"]

    if 'guspec' in df.columns:
        result = check_against_ldms_with_guspec(df)
    elif 'guspec_core' in df.columns:
        result = check_against_ldms_with_guspec_core(df)
    else:
        result = "guspec/guspec_core col missing"
    return result

def check_against_ldms_with_guspec_core(df):
    ldms = pd.DataFrame()
    for i, row in df[['network','protocol']].drop_duplicates().iterrows():
        s = get_ldms(row.network, row.protocol)
        ldms = pd.concat([ldms, s])
    ldms = ldms.rename(columns=constants.LDMS_RELABEL_DICT)
    ldms["drawdt"] = ldms.apply(
        lambda x: datetime.date(x.drawdy, x.drawdm, x.drawdd).isoformat(), axis=1
    )
    ldms = ldms.drop(columns=["drawdy", "drawdm", "drawdd"])
    ldms['guspec_core'] = ldms.guspec.str.rpartition("-")[0]
    ldms = ldms.loc[ldms.guspec_core.isin(df.guspec_core)]

    ldms.protocol = ldms.protocol.astype(int)
    ldms.ptid = ldms.ptid.astype(int)
    ldms.visitno = ldms.visitno.astype(float)
    usecols = list(set(ldms.columns).intersection(df.columns))

    ldms = ldms[usecols]
    ldms = ldms.drop_duplicates()

    # cast df types
    df.ptid = df.ptid.astype(int)
    df.visitno = df.visitno.astype(float)
    df.protocol = df.protocol.astype(int)
    df.drawdt = pd.to_datetime(df['drawdt']).astype(str)
    for col in ['guspec_core','spec_primary','spec_additive','spec_derivative']:
        if col in df.columns:
            df[col] = df[col].astype(str)
            ldms[col] = ldms[col].astype(str)

    df_ldms = df[usecols].drop_duplicates()

    check_outer = df_ldms.merge(ldms, on=usecols, how='outer')
    check_inner = df_ldms.merge(ldms, on=usecols, how='inner')
    # nancount = check_outer.sdmc_processing_datetime.isna().sum() + check_inner.sdmc_processing_datetime.isna().sum()
    # print(len(df_ldms), len(check_outer), len(check_inner))
    if len(df_ldms)!=len(check_outer) or len(df_ldms)!=len(check_inner):
        return "mismatch"
    else:
        df_ldms = df_ldms.sort_values(by=usecols).reset_index(drop=True)
        ldms = ldms.sort_values(by=usecols).reset_index(drop=True)

        equality_check = df_ldms.compare(ldms)
        if len(equality_check) == 0:
            return "looks good"
        else:
            return "mismatch"

def check_against_ldms_with_guspec(df):
    ldms = pd.DataFrame()
    for i, row in df[['network','protocol']].drop_duplicates().iterrows():
        s = get_ldms(row.network, row.protocol)
        ldms = pd.concat([ldms, s])
    ldms = ldms.rename(columns=constants.LDMS_RELABEL_DICT)
    ldms["drawdt"] = ldms.apply(
        lambda x: datetime.date(x.drawdy, x.drawdm, x.drawdd).isoformat(), axis=1
    )
    ldms = ldms.drop(columns=["drawdy", "drawdm", "drawdd"])
    ldms = ldms.loc[ldms.guspec.isin(df.guspec)]

    ldms.protocol = ldms.protocol.astype(int)
    ldms.ptid = ldms.ptid.astype(int)
    ldms.visitno = ldms.visitno.astype(float)

    usecols = list(set(ldms.columns).intersection(df.columns))
    ldms = ldms[usecols]
    ldms = ldms.drop_duplicates()

    # cast df types
    df.ptid = df.ptid.astype(int)
    df.visitno = df.visitno.astype(float)
    df.protocol = df.protocol.astype(int)
    df.drawdt = pd.to_datetime(df['drawdt']).astype(str)
    for col in ['guspec','spec_primary','spec_additive','spec_derivative']:
        if col in df.columns:
            df[col] = df[col].astype(str)
            ldms[col] = ldms[col].astype(str)

    df_ldms = df[usecols].drop_duplicates()

    check_outer = df_ldms.merge(ldms, on=usecols, how='outer')
    check_inner = df_ldms.merge(ldms, on=usecols, how='inner')
    # nancount = check_outer.sdmc_processing_datetime.isna().sum() + check_inner.sdmc_processing_datetime.isna().sum()
    # print(len(df_ldms), len(check_outer), len(check_inner))
    if len(df_ldms)!=len(check_outer) or len(df_ldms)!=len(check_inner):
        return "mismatch"
    else:
        df_ldms = df_ldms.sort_values(by=usecols).reset_index(drop=True)
        ldms = ldms.sort_values(by=usecols).reset_index(drop=True)

        equality_check = df_ldms.compare(ldms)
        if len(equality_check) == 0:
            return "looks good"
        else:
            return "mismatch"

def get_output_path(yamlpath):
    y = read_yaml(yamlpath)
    savedir = y['savedir']
    files = os.listdir(savedir)
    outputpath = np.sort([i for i in files if 'process' in i.lower() and '.txt' in i.lower()])[-1]
    if savedir[-1]!="/":
        savedir += "/"
    return savedir + outputpath

def find_endpoints(d: str, l: list) -> list:
    """
    for everything in dir d,
    return if it's an endpoint (non-dir),
    in addition to the contents of l, as a list.
    """
    if os.path.isdir(d):
        for c in os.listdir(d):
            l = find_endpoints(d + "/" + c, l)
        return l
    else:
        return l + [d]

def read_yaml(yaml_path: str, add_path:bool = False) -> dict:
    """
    given yaml path
    optionally add yaml_path to yaml dict
    return yaml dict
    """
    with open(yaml_path, 'r') as file:
        yaml_dict = yaml.safe_load(file)
    if yaml_dict is None:
        # print(f"{yaml_path} is currently empty; please fill in")
        return {"yaml_path": yaml_path}
    else:
        yaml_dict["yaml_path"] = yaml_path
        return yaml_dict

def get_output_path_from_yaml(yaml_dict: dict) -> list:
    files = os.listdir(yaml_dict['savedir'])
    if not isinstance(yaml_dict['output_prefix'], list):
        yaml_dict['output_prefix'] = [yaml_dict['output_prefix']]
    fnames = []
    for prefix in yaml_dict['output_prefix']:
        fnames += [np.sort([f for f in files if prefix in f])[-1]]
    output_path = yaml_dict['savedir']
    if output_path[-1] != "/":
        output_path += "/"
    return [output_path + f for f in fnames]




if __name__=="__main__":
    main()
