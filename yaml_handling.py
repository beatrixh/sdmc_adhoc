import yaml
import pandas as pd
import numpy as np
import os

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

def overwrite_yaml_protocols(network: str, protocols: list) -> None:
    """
    given a network and the corresponding protocols,
    overwrite in sdmc-adhoc/constants.yaml the list of
    all protocols for the specified network
    """
    protocols = np.array(protocols)
    protocols = protocols[~np.isnan(protocols)].astype(int)
    protocols = np.unique(protocols).tolist()

    yamlpath = os.path.dirname(__file__) + "/constants.yaml"
    with open(yamlpath, 'r') as file:
        yamldict = yaml.safe_load(file)

    yamldict[f"{network.lower()}_protocols"] = protocols

    with open(yamlpath, 'w') as outfile:
        yaml.dump(yamldict, outfile, default_flow_style=False, sort_keys=False)


def add_to_yaml_protocols(network: str, new_protocols: list) -> None:
    """
    given a network and corresponding protocols,
    add to the list of protocols in sdmc-adhoc/constants.yaml
    if the new protocols aren't already tracked
    """
    # read in yaml
    yamlpath = os.path.dirname(__file__) + "/constants.yaml"
    with open(yamlpath, 'r') as file:
        yamldict = yaml.safe_load(file)

    # add new protocols to list
    protocols = yamldict[f"{network.lower()}_protocols"]
    protocols += new_protocols

    # make sure list is unique, contains no nans
    protocols = np.array(protocols)
    protocols = protocols[~np.isnan(protocols)].astype(int)
    protocols = np.unique(protocols).tolist()

    # add updated list to yaml dict
    yamldict[f"{network.lower()}_protocols"] = protocols

    # save to constants.yaml
    with open(yamlpath, 'w') as outfile:
        yaml.dump(yamldict, outfile, default_flow_style=False, sort_keys=False)

def get_protocol_from_output(datapath: str) -> dict:
    """
    given a path to a datafile of outputs
    that contains "network" and "protocol" columns
    (no specific capitalization),
    return a dictionary with keys hvtn/covpn
    corresponding to a list[int] of protocols
    """
    if not os.path.exists(datapath):
        raise Exception(f"This is a dead path: {datapath}")
    if datapath[-3:]=="txt":
        data = pd.read_csv(datapath, sep="\t")
    elif datapath[-3:]=="csv":
        data = pd.read_csv(datapath)
    else:
        raise Exception(f"Filepath neither txt nor csv: {datapath}")

    #convert to lowercase just in case
    data.columns = [i.lower() for i in data.columns]
    hvtn = data.loc[
        (data.network.str.lower()=="hvtn") & (data.protocol.notna())
        ].protocol.astype(int).unique().tolist()
    covpn = data.loc[
        (data.network.str.lower()=="covpn") & (data.protocol.notna())
        ].protocol.astype(int).unique().tolist()

    return {"hvtn": hvtn, "covpn": covpn}

def get_protocols_from_old_outputs() -> dict:
    """
    get a list of all applicable protocols for each network
    in old (pre-2024) outputs, and return as a dict
    with keys "hvtn"/"covpn"
    """
    hvtn, covpn = set(), set()
    for path in OUTPUT_REGISTRY_OLD:
        d = get_protocol_from_output(path)
        hvtn.update(d["hvtn"])
        covpn.update(d["covpn"])
    return {"hvtn": list(hvtn), "covpn": list(covpn)}


def get_protocols_from_new_outputs() -> dict:
    """
    get a list of all applicable protocols for each network
    in new (2024+) outputs, and return as a dict
    with keys "hvtn"/"covpn"
    """
    hvtn, covpn = set(), set()
    yamls = [
        f for f in find_endpoints('/home/bhaddock/repos/sdmc-adhoc/processing_scripts', l=[]) if f.split("/")[-1]=="paths.yaml"
    ]
    for y in yamls:
        yaml_dict = read_yaml(y, add_path=False)
        output_filepaths = get_output_path_from_yaml(yaml_dict)
        for path in output_filepaths:
            d = get_protocol_from_output(path)
            hvtn.update(d["hvtn"])
            covpn.update(d["covpn"])
    return {"hvtn": list(hvtn), "covpn": list(covpn)}

if __name__=='__main__':
    yamlpath = os.path.dirname(__file__) + "/constants.yaml"
    yamldict = read_yaml(yamlpath)
    OUTPUT_REGISTRY_OLD = yamldict["OUTPUT_REGISTRY_OLD"]

    # old_protocols = get_protocols_from_old_outputs()
    new_protocols = get_protocols_from_new_outputs()

    # hvtn = set(old_protocols["hvtn"]).union(new_protocols["hvtn"])
    # covpn = set(old_protocols["covpn"]).union(new_protocols["covpn"])

    print(f"found hvtn protocols: {new_protocols['hvtn']}")
    print(f"also covpn protocols: {new_protocols['covpn']}")
    # overwrite_yaml_protocols("hvtn", list(hvtn))
    # overwrite_yaml_protocols("covpn", list(covpn))
