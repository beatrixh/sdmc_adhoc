## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 01/17/2025
# Purpose: Pull a DataFrame each column name that's been in ad hoc processing,
# and count of its frequency
## ---------------------------------------------------------------------------##
import yaml
import pandas as pd
import numpy as np
import os
import sys

SAVEPATH = sys.argv[1]
REMOVE_THESE = [
    '/home/bhaddock/repos/sdmc-adhoc/processing_scripts/TEMPLATE/LAB_ASSAY/paths.yaml',
    '/home/bhaddock/repos/sdmc-adhoc/processing_scripts/TEMPLATE/TEST_ASSAY/paths.yaml',
]

def main():
    # find all files in processing_scripts dir
    endpoints = find_endpoints('/home/bhaddock/repos/sdmc-adhoc/processing_scripts', l=[])

    # subset to only .yamls, then read them in
    yamls = [i for i in endpoints if i[-4:]=='yaml']

    # hand-pick ones to exclude
    yamls = list(set(yamls).difference(REMOVE_THESE))

    # read in yamls
    yamls = [read_yaml(y, add_path=False) for y in yamls]

    # pull output paths from yamls
    output_paths = [get_output_path_from_yaml(y) for y in yamls]

    # filter out output paths that broke
    non_matches = [i for i in output_paths if 'NO MATCH' in i]
    if len(non_matches) > 0:
        print(f"The following output paths weren't located; go fix yamls: {non_matches}")
    output_paths = [i for i in output_paths if 'NO MATCH' not in i]

    # flatten list
    output_paths = np.concatenate(output_paths).tolist()

    # pull list of all column names used
    colnames = []
    for o in output_paths:
        df = pd.read_csv(o, sep="\t")
        colnames += list(df.columns)

    # get freq of colnames
    colname_counts = pd.Series(colnames).value_counts().reset_index()
    colname_counts.columns = ['col', 'N']

    colname_counts.to_csv(SAVEPATH, index=False)
    print(colname_counts)

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
    try:
        if not isinstance(yaml_dict['output_prefix'], list):
            yaml_dict['output_prefix'] = [yaml_dict['output_prefix']]
        fnames = []
        for prefix in yaml_dict['output_prefix']:
            fnames += [np.sort([f for f in files if prefix.lower() in f.lower() and 'process' in f.lower() and 'notes' not in f.lower()])[-1]]
        output_path = yaml_dict['savedir']
        if output_path[-1] != "/":
            output_path += "/"
        return [output_path + f for f in fnames]
    except:
        try:
            return yaml_dict['savedir'] + 'NO MATCH'
        except:
            return 'NO MATCH'

if __name__=="__main__":
    main()
