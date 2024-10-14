import pandas as pd, numpy as np
import os
import yaml
import datetime
from typing import List

import sdmc_tools.constants as constants
from yaml_handling import find_endpoints, read_yaml, get_output_path_from_yaml

## check if ldms has been updated today ------------------------------------- ##
def was_ldms_updated(NETWORK: str) -> bool:
    """
    given a network
    - if ldms has been updated since the last weekday, return True
    - if ldms hasn't been updated since last weekday, return False
    """
    # get today's timestamp
    today = datetime.date.today()
    timestamp = today.strftime("%Y%m%d")
    if today.weekday() == 0: # if it's a monday
        last_weekday = today - datetime.timedelta(days = 3)
    else: #any other day. note this is assuming the script isnt run sat or sun
        last_weekday = today - datetime.timedelta(days = 1)
    # get timestamp that LDMS most recently updated
    if NETWORK=="HVTN":
        LAST_TOUCHED = os.stat(constants.LDMS_PATH_HVTN).st_mtime
    elif NETWORK=="CoVPN":
        LAST_TOUCHED = os.stat(constants.LDMS_PATH_COVPN).st_mtime
    else:
        print("NETWORK MUST BE EITHER HVTN OR CoVPN")

    LAST_TOUCHED = datetime.datetime.fromtimestamp(LAST_TOUCHED)
    return LAST_TOUCHED.date() > last_weekday

## save todays ldms --------------------------------------------------------- ##
def save_todays_ldms(PROTOCOLS: List[str], NETWORK: str) -> None:
    """
    INPUTS
     - PROTOCOLS: list of PROTOCOLS (numeric) from one network
     - NETWORK: string ("HVTN" or "CoVPN") corresponding to above PROTOCOLS
    """
    if NETWORK=="HVTN":
        ldms = pd.read_csv(constants.LDMS_PATH_HVTN,
                           usecols=constants.STANDARD_COLS,
                           dtype=constants.LDMS_DTYPE_MAP
                           )
    elif NETWORK=="CoVPN":
        ldms = pd.read_csv(constants.LDMS_PATH_COVPN,
                           usecols=constants.STANDARD_COLS,
                           dtype=constants.LDMS_DTYPE_MAP
                           )
    else:
        print(f"{NETWORK} IS AN INVALID NETWORK SELECTION")
    ldms = ldms.loc[ldms.lstudy.isin(PROTOCOLS)]
    timestamp = datetime.date.today().strftime("%Y%m%d")

    network_path = f"/networks/vtn/lab/SDMC_labscience/studies/{NETWORK}/"
    for protocol in PROTOCOLS:
        protocol_path = network_path + PROTOCOL_DIRNAME_MAP[NETWORK][int(protocol)]
        if not os.path.exists(protocol_path):
            print(f"creating protocol dir for {int(protocol)}")
            os.mkdir(protocol_path)
        specimen_path = protocol_path + "/specimens"
        if not os.path.exists(specimen_path):
            print(f"creating specimens dir for {int(protocol)}")
            os.mkdir(specimen_path)
        ldms_feed_path = specimen_path + "/ldms_feed"
        if not os.path.exists(ldms_feed_path):
            print(f"creating ldms_feed dir for {int(protocol)}")
            print(ldms_feed_path)
            os.mkdir(ldms_feed_path)

        savepath = ldms_feed_path + f"/{NETWORK.lower()}.ldms{int(protocol)}.{timestamp}.csv"
        if not os.path.exists(savepath):
            # print(f"CREATING THE FOLLOWING: {savepath}")
            ldms.loc[ldms.lstudy==protocol].to_csv(savepath, index=False)
        else:
            print(f"{savepath} already exists")

## delete old ldms ---------------------------------------------------------- ##
def delete_old_ldms(PROTOCOL: str, NETWORK: str) -> None:
    PROTOCOL = int(PROTOCOL)
    N = len(NETWORK) + len(".ldms") + len(str(PROTOCOL))
    feed_dir = f"/networks/vtn/lab/SDMC_labscience/studies/{NETWORK}/{PROTOCOL_DIRNAME_MAP[NETWORK][int(PROTOCOL)]}/specimens/ldms_feed/"
    things_in_dir = os.listdir(feed_dir)
    # sort from oldest to most recent, and only take the ones that look like hvtn.ldmsPROTOCOL
    applicable_things_in_dir = np.sort([i for i in things_in_dir if i[:N] == f'{NETWORK.lower()}.ldms{PROTOCOL}'])

    # delete all but the two most recent
    for file in applicable_things_in_dir[:-2]:
        # print(f"DELETING THE FOLLOWING: {file}\n")
        os.remove(feed_dir + file)

## detect_ldms_diffs -------------------------------------------------------- ##
def detect_ldms_diffs(PROTOCOL: str, NETWORK: str) -> None:
    """
    INPUT
        - protocol number (int)
        - network (HVTN or COVPN) associated with protocol
    FUNCTION: determines if yesterday / todays ldms has updates in any of the columns we're interested in
    """
    # print(f"\n{NETWORK}{PROTOCOL}: Reading in data")
    old, new = get_ldms_subset(PROTOCOL, NETWORK)

    # if either of these are empty, we have nothing to compare ----------------#
    if old.empty or new.empty:
        return

    # the columns have to match -- otherwise have big problems ----------------#
    # print(f"{NETWORK}{PROTOCOL}: Checking for col diff")
    COL_DIFF = set(old.columns).symmetric_difference(new.columns)
    if len(COL_DIFF) > 0:
        print(f"{NETWORK}{PROTOCOL}: COLUMNS DON'T MATCH. THE FOLLOWING NOT SHARED: {COL_DIFF}")
        return

    # do we have the same set of guspecs --------------------------------------#
    guspecs_removed = set(old.guspec).difference(new.guspec)
    if guspecs_removed:
        print(f"{NETWORK}{PROTOCOL}: THE FOLLOWING GUSPECS HAVE BEEN REMOVED FROM LDMS: {guspecs_removed}")
        handle_affected_jobs(guspecs_removed)
    guspecs_added = set(new.guspec).difference(old.guspec)
    if guspecs_added:
        print(f"THE FOLLOWING GUSPECS HAVE BEEN ADDED TO LDMS {guspecs_added}")
        handle_affected_jobs(guspecs_added)
    # if not guspecs_removed and not guspecs_added:
        # print(f"{NETWORK}{PROTOCOL}: No guspecs added or removed")

    # investigate all rows with shared guspec----------------------------------#
    # print(f"{NETWORK}{PROTOCOL}: Comparing all rows with shared guspecs")
    shared_guspecs = set(old.guspec).intersection(new.guspec)

    new = new.loc[new.guspec.isin(shared_guspecs)].drop_duplicates()
    old = old.loc[old.guspec.isin(shared_guspecs)].drop_duplicates()

    new_guspec_count = new.guspec.value_counts().to_frame()
    old_guspec_count = old.guspec.value_counts().to_frame()

    comp = new_guspec_count.merge(
        old_guspec_count,
        how='outer',
        left_index=True,
        right_index=True
    )

    GUSPEC_DIFFS = list(comp.loc[comp.count_x!=comp.count_y].index)

    # report on any guspecs corresponding to different row counts -------------#
    if GUSPEC_DIFFS:
        print(f"{NETWORK}{PROTOCOL}: DIFFERENT NUMBER OF ROWS PER SHARED GUSPEC")
        print(f"DIFFERENCES WITH THE FOLLOWING GUSPECS: {GUSPEC_DIFFS}")
        handle_affected_jobs(GUSPEC_DIFFS)

        print("Comparing matching rows:")
        GUSPEC_MATCHES = list(comp.loc[comp.count_x==comp.count_y].index)
        new = new.loc[new.guspec.isin(GUSPEC_MATCHES)]
        old = old.loc[old.guspec.isin(GUSPEC_MATCHES)]
    # else:
    #     print("For shared guspecs, we have the exact count of each guspec")

    # report on guspecs with matching row counts ------------------------------#
    # print("comparing rows with shared guspec / guspec count:")
    # print(f"{NETWORK}{PROTOCOL}: Sorting dataframes")

    old = old.sort_values(
        by=['guspec','txtpid', 'drawdm', 'drawdd', 'drawdy', 'vidval',
        'lstudy', 'primstr', 'addstr', 'dervstr'],
        ignore_index=True
    )
    new = new.sort_values(
        by=['guspec','txtpid', 'drawdm', 'drawdd', 'drawdy', 'vidval',
        'lstudy', 'primstr', 'addstr', 'dervstr'],
        ignore_index=True
    )

    # catch any rows that have mismatch
    # print(f"{NETWORK}{PROTOCOL}: Storing dataframe diff")
    diff = new.compare(old)
    if len(diff) > 0:
        print(f"CHANGES DETECTED IN SHARED GUSPECS")
        AFFECTED_GUSPECS = list(new.iloc[diff.index].guspec.unique())
        print(f"AFFECTED GUSPECS: {AFFECTED_GUSPECS}")
        handle_affected_jobs(AFFECTED_GUSPECS)
    # else:
    #     print(f"NO DIFF DETECTED FOR {NETWORK}{PROTOCOL}")


# helpers for detect_ldms_diffs ----------------------------------------------##
def get_ldms_subset(PROTOCOL: str, NETWORK: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Given a protocol and network
    return the prior saved ldms [CURRENTLY DEFINED ONLY AS LAST WEEKDAY] and today's ldms
    """
    today = datetime.date.today()
    timestamp = today.strftime("%Y%m%d")

    feed_dir = f"/networks/vtn/lab/SDMC_labscience/studies/{NETWORK}/{PROTOCOL_DIRNAME_MAP[NETWORK][int(PROTOCOL)]}/specimens/ldms_feed/"
    fname = f"{NETWORK.lower()}.ldms{PROTOCOL}.{timestamp}.csv"

    if os.path.exists(feed_dir + fname):
        new = pd.read_csv(feed_dir + fname, dtype=constants.LDMS_DTYPE_MAP)
        files = os.listdir(feed_dir)
        if len(files) > 1:
            prev_fname = np.sort(os.listdir(feed_dir))[-2]
            old = pd.read_csv(feed_dir + prev_fname, dtype=constants.LDMS_DTYPE_MAP)
            # print(f"{NETWORK}{PROTOCOL}: Comparing {fname} and {prev_fname}")
        else:
            print(f"NO PRIOR LDMS SAVED FOR {NETWORK}{PROTOCOL}")
            old = pd.DataFrame()
    else:
        print(f"NO LDMS SAVED FOR {NETWORK}{PROTOCOL} TODAY")
        new = pd.DataFrame()

    return old, new

def handle_affected_jobs(guspecs: list) -> None:
    handle_affected_jobs_old(guspecs)
    handle_affected_jobs_new(guspecs)

def handle_affected_jobs_old(guspecs: list) -> None: #report_outputs_corresponding_to_guspecs
    """
    given a list of guspecs with changes in ldms,
    report on any that are used by older (non-auto rerunnable) jobs
    """
    affected_outputs = set()
    for guspec in guspecs:
        if guspec in GUSPEC_TO_OUTPUT_PATH_OLD.keys():
            affected_outputs.update(GUSPEC_TO_OUTPUT_PATH_OLD[guspec])
    if len(affected_outputs) > 0:
        print(f"AFFECTED OUTPUTS: {list(AFFECTED_OUTPUTS)}")
    else:
        print("\nNO PRE-2024 OUTPUTS AFFECTED")

def handle_affected_jobs_new(guspecs: list) -> None:
    """
    given a list guspecs with changes in ldms,
    return a list of yaml dicts corresponding to jobs that used that guspec
    """
    # pull all yaml_dicts associated with jobs
    print("Checking to see if new (2024+) outputs affected\n")
    # print("Finding yamls")
    yamls = [
        f for f in find_endpoints('/home/bhaddock/repos/sdmc-adhoc/processing_scripts', l=[]) if f.split("/")[-1]=="paths.yaml"
    ]
    # print("yamls obtained")
    yamls = [read_yaml(p, add_path=True) for p in yamls]
    affected_job_yamls = []
    # print("looping through guspecs, checking if in yamls")
    # print(f"guspecs to loop through: {guspecs}\n")
    # print(f"yamls to loop through: {yamls}")
    for guspec in guspecs:
        for y in yamls:
            if "guspecs" not in y.keys():
                print("adding guspecs to yaml\n")
                save_guspecs_to_yaml(y["yaml_path"])
                y = read_yaml(y["yaml_path"], add_path=True)
            if guspec in y["guspecs"]:
                print(f"found guspec in an output: {guspec}")
                affected_job_yamls += [y]
    if len(affected_job_yamls) > 0:
        print("\nNEW OUTPUTS AFFECTED")
        for y in affected_job_yamls:
            output_path = get_output_path_from_yaml(y)
            print(f"The following output file affected: {output_path}")
            rerun_affected_job(y)
    else:
        print("\nNO NEW (2024+) OUTPUTS AFFECTED")

def save_guspecs_to_yaml(yaml_path: str) -> None:
    """
    given a filepath to a yaml that contains
    - a savedir
    - an output prefix(es)
    save the guspecs from the output to the yaml under key "guspecs"
    """
    # print("reading in yaml\n")
    yaml_dict = read_yaml(yaml_path, add_path=False)
    # print("getting output path from yaml")
    output_filepaths = get_output_path_from_yaml(yaml_dict)
    guspecs = []
    if not isinstance(output_filepaths, list):
        output_filepaths = [output_filepaths]
    for path in output_filepaths:
        if path[-3:]=="txt":
            data = pd.read_csv(path, sep="\t")
        elif path[-3:]=="csv":
            data = pd.read_csv(path)
        else:
            raise Exception(f"PUT EXCEPTION HERE; FILEPATH {path} WRONG")
        if 'guspec' in data.columns:
            guspecs += data.guspec.tolist()
    yaml_dict['guspecs'] = list(set(guspecs))
    with open(yaml_path, 'w') as outfile:
        yaml.dump(yaml_dict, outfile, default_flow_style=False, sort_keys=False)

def rerun_affected_job(affected_job_yaml: str) -> None:
    """
    given a yaml, report that we're trying to regenerate
    the corresponding output and rerun the corresponding script
    """
    yaml_path = affected_job_yaml["yaml_path"]
    # yaml = '/home/bhaddock/repos/sdmc-adhoc/processing_scripts/HVTN128/PKSMC/paths.yaml'
    script_path = yaml_path[len("/home/bhaddock/repos/sdmc-adhoc/"):-len("/paths.yaml")] + "/process_data.py"
    print(f"Rerunning {script_path}")
    package = script_path[:-3].replace("/",".")
    try:
        rerun = getattr(__import__(package, fromlist=["main"]), "main")
        rerun()
    except:
        print(f"ERROR trying to run {script_path}")

## pull in constants -------------------------------------------------------- ##
yamlpath = os.path.dirname(__file__) + "/constants.yaml"
yamldict = read_yaml(yamlpath)

PROTOCOL_DIRNAME_MAP = yamldict["PROTOCOL_DIRNAME_MAP"]
GUSPEC_TO_OUTPUT_PATH_OLD = yamldict["GUSPEC_TO_OUTPUT_PATH_OLD"]


# def read_in_output(path):
#     """
#     INPUT: filepath to csv or txt of data
#     OUTPUT: data contained at filepath, columns standardized
#     """
#     filetype = path.rpartition(".")[-1]
#     if filetype =="txt":
#         data = pd.read_csv(path, sep="\t")
#     elif filetype == "csv":
#         data = pd.read_csv(path)
#     else:
#         return f"UNRECOGNIZED FILETYPE: {filetype}"
#     # if we were able to read in the data
#     data.columns = [i.lower() for i in data.columns]
#     return data

# def get_output_filepath_from_yaml(yamlpath: str) -> list:
#     with open(yamlpath, 'r') as file:
#         yaml_dict = yaml.safe_load(file)
#     output_dir = yaml_dict["savedir"]
#     if output_dir[-1]!="/":
#         output_dir += "/"
#     #TODO: add try/excepts for if there's a nonmatch or this is missing
#     try:
#         prefixes = yaml_dict["output_prefix"]
#         if not isinstance(prefixes, list):
#             prefixes = [prefixes]
#         fnames = []
#         for prefix in prefixes:
#             matches = [file for file in os.listdir(output_dir) if file[:len(prefix)]==prefix]
#             fnames += [output_dir + np.sort(matches)[-1]]
#         return fnames
#     except:
#         return [output_dir + "MISSING"]
