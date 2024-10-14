import pandas as pd
import numpy as np
import datetime
import yaml
import os
import sdmc_tools.constants as constants

yamlpath = '/home/bhaddock/repos/sdmc-adhoc/constants.yaml'
with open(yamlpath, 'r') as file:
    yamldict = yaml.safe_load(file)

def main():
    # this is only checking for changes in new outputs
    endpoints = find_endpoints("/home/bhaddock/repos/sdmc-adhoc/processing_scripts", l=[])
    yams = [i for i in endpoints if i[-4:]=="yaml"]
    output_paths = [get_output_path(i) for i in yams]

    for o in output_paths:
        try:
            check_for_changes(o)
        except:
            print(f"Encountered issue with {o}")
def get_output_path(yamlpath):
    y = read_yaml(yamlpath)
    savedir = y['savedir']
    files = os.listdir(savedir)
    outputpath = np.sort([i for i in files if 'process' in i.lower() and '.txt' in i.lower()])[-1]
    if savedir[-1]!="/":
        savedir += "/"
    return savedir + outputpath

def get_ldms(network, protocol):
    dirname = yamldict['PROTOCOL_DIRNAME_MAP'][network][int(protocol)]
    feeddir = f"/networks/vtn/lab/SDMC_labscience/studies/{network}/{dirname}/specimens/ldms_feed/"
    fname = np.sort(os.listdir(feeddir))[-1]
    ldms_path = feeddir + fname
    return pd.read_csv(ldms_path, dtype=constants.LDMS_DTYPE_MAP)

def check_for_changes(o):
    df = pd.read_csv(o, sep="\t")
    if 'guspec' in [i.lower() for i in df.columns]:
        network = df.network.unique()[0]
        protocol = df.protocol.unique()[0]

        ldms = get_ldms(network, protocol)
        ldms = ldms.rename(columns=constants.LDMS_RELABEL_DICT)
        ldms["drawdt"] = ldms.apply(
            lambda x: datetime.date(x.drawdy, x.drawdm, x.drawdd).isoformat(), axis=1
        )
        ldms = ldms.drop(columns=["drawdy", "drawdm", "drawdd"])
        ldms = ldms.loc[ldms.guspec.isin(df.guspec)]

        ldms.protocol = ldms.protocol.astype(int)
        ldms.ptid = ldms.ptid.astype(int)
        ldms = ldms.drop_duplicates()

        check = df.merge(ldms, on=ldms.columns.tolist(), how='right')
        if len(check)!=len(df):
            print(f"SOMETHING NOT MATCHING, {o}")
            print(f"len(check): {len(check)}; len(df): {len(df)}")

        nancount = check.sdmc_processing_datetime.isna().sum()
        if nancount > 0:
            print(f"SOMETHING NOT MATCHING, {o}")
            print(f"sdmc_processing_datetime nan count: {nancount}")
    else:
        print(f"no guspec col in {o}")

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
        print(f"{yaml_path} is currently empty; please fill in")
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
