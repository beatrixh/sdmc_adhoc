import pandas as pd, numpy as np
import os
import sdmc_tools.constants as constants
from ldms_monitoring_constants import *
import datetime
from typing import List

def was_ldms_updated(NETWORK):
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

def save_todays_ldms(PROTOCOLS: List[str], NETWORK: str):
    """
    INPUTS
     - PROTOCOLS: list of PROTOCOLS (numeric) from one network
     - NETWORK: string ("HVTN" or "CoVPN") corresponding to above PROTOCOLS
    """
    if NETWORK=="HVTN":
        ldms = pd.read_csv(constants.LDMS_PATH_HVTN,
                           usecols=constants.STANDARD_COLS,
                           dtype=dtype_map
                           )
    elif NETWORK=="CoVPN":
        ldms = pd.read_csv(constants.LDMS_PATH_COVPN,
                           usecols=constants.STANDARD_COLS,
                           dtype=dtype_map
                           )
    else:
        print(f"{NETWORK} IS AN INVALID NETWORK SELECTION\n")
    ldms = ldms.loc[ldms.lstudy.isin(PROTOCOLS)]
    timestamp = datetime.date.today().strftime("%Y%m%d")

    network_path = f"/networks/vtn/lab/SDMC_labscience/studies/{NETWORK}/"
    for protocol in PROTOCOLS:
        protocol_path = network_path + PROTOCOL_DIRNAME_MAP[NETWORK][int(protocol)]
        if not os.path.exists(protocol_path):
            print(f"creating protocol dir for {int(protocol)}\n")
            os.mkdir(protocol_path)
        specimen_path = protocol_path + "/specimens"
        if not os.path.exists(specimen_path):
            print(f"creating specimens dir for {int(protocol)}\n")
            os.mkdir(specimen_path)
        ldms_feed_path = specimen_path + "/ldms_feed"
        if not os.path.exists(ldms_feed_path):
            print(f"creating ldms_feed dir for {int(protocol)}\n")
            print(ldms_feed_path)
            os.mkdir(ldms_feed_path)

        savepath = ldms_feed_path + f"/{NETWORK.lower()}.ldms{int(protocol)}.{timestamp}.csv"
        if not os.path.exists(savepath):
            print(f"CREATING THE FOLLOWING: {savepath}\n")
            ldms.loc[ldms.lstudy==protocol].to_csv(savepath, index=False)
        else:
            print(f"{savepath} already exists\n")

def delete_old_ldms(PROTOCOL, NETWORK):
    PROTOCOL = int(PROTOCOL)
    N = len(NETWORK) + len(".ldms") + len(str(PROTOCOL))
    feed_dir = f"/networks/vtn/lab/SDMC_labscience/studies/{NETWORK}/{PROTOCOL_DIRNAME_MAP[NETWORK][int(PROTOCOL)]}/specimens/ldms_feed/"
    things_in_dir = os.listdir(feed_dir)
    # sort from oldest to most recent, and only take the ones that look like hvtn.ldmsPROTOCOL
    applicable_things_in_dir = np.sort([i for i in things_in_dir if i[:N] == f'{NETWORK.lower()}.ldms{PROTOCOL}'])

    # delete all but the two most recent
    for file in applicable_things_in_dir[:-2]:
        print(f"DELETING THE FOLLOWING: {file}\n")
        os.remove(feed_dir + file)

def detect_ldms_diffs(PROTOCOL, NETWORK):
    """
    INPUT
        - protocol number (int)
        - network (HVTN or COVPN) associated with protocol
    FUNCTION: determines if yesterday / todays ldms has updates in any of the columns we're interested in
    """
    print(f"{NETWORK}{PROTOCOL}: reading in data")
    old, new = get_ldms_subset(PROTOCOL, NETWORK)

    # if either of these are empty, we have nothing to compare ----------------#
    if old.empty or new.empty:
        return

    # the columns have to match -- otherwise have big problems ----------------#
    print(f"{NETWORK}{PROTOCOL}: checking for col diff")
    COL_DIFF = set(old.columns).symmetric_difference(new.columns)
    if len(COL_DIFF) > 0:
        print(f"COLUMNS DON'T MATCH. THE FOLLOWING NOT SHARED: {COL_DIFF}\n")
        return

    # do we have the same set of guspecs --------------------------------------#
    guspecs_removed = set(old.guspec).difference(new.guspec)
    if guspecs_removed:
        print(f"THE FOLLOWING GUSPECS HAVE BEEN REMOVED FROM LDMS: {guspecs_removed}\n")
        handle_affected_jobs(guspecs_removed)
    guspecs_added = set(new.guspec).difference(old.guspec)
    if guspecs_added:
        print(f"THE FOLLOWING GUSPECS HAVE BEEN ADDED TO LDMS {guspecs_added}\n")
        handle_affected_jobs(guspecs_added)
    if not guspecs_removed and not guspecs_added:
        print(f"{NETWORK}{PROTOCOL}: No guspecs added or removed")

    # investigate all rows with shared guspec----------------------------------#
    print(f"{NETWORK}{PROTOCOL}: Comparing all rows with shared guspecs")
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

        print("comparing matching rows:")
        GUSPEC_MATCHES = list(comp.loc[comp.count_x==comp.count_y].index)
        new = new.loc[new.guspec.isin(GUSPEC_MATCHES)]
        old = old.loc[old.guspec.isin(GUSPEC_MATCHES)]
    else:
        print("\nfor shared guspecs, we have the exact count of each guspec")

    # report on guspecs with matching row counts ------------------------------#
    print("comparing rows with shared guspec / guspec count:")
    print(f"{NETWORK}{PROTOCOL}: sorting dataframes")

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
    print(f"{NETWORK}{PROTOCOL}: storing dataframe diff")
    diff = new.compare(old)
    if len(diff) > 0:
        print(f"CHANGES DETECTED IN SHARED GUSPECS\n")
        AFFECTED_GUSPECS = list(new.iloc[diff.index].guspec.unique())
        print(f"AFFECTED GUSPECS: {AFFECTED_GUSPECS}\n")
        handle_affected_jobs(AFFECTED_GUSPECS)
    else:
        print(f"NO DIFF DETECTED FOR {NETWORK}{PROTOCOL}\n")

# helpers for detect_ldms_diffs ----------------------------------------------##
def get_ldms_subset(PROTOCOL, NETWORK):
    """
    Given a protocol and network
    return the prior saved ldms [CURRENTLY DEFINED ONLY AS LAST WEEKDAY] and today's ldms
    """
    today = datetime.date.today()
    timestamp = today.strftime("%Y%m%d")
    # if today.weekday() == 0: # if it's a monday
    #     yesterday = today - datetime.timedelta(days = 3)
    # else: #any other day. note this is assuming the script isnt run sat or sun
    #     yesterday = today - datetime.timedelta(days = 1)
    # prev_timestamp = yesterday.strftime("%Y%m%d")

    feed_dir = f"/networks/vtn/lab/SDMC_labscience/studies/{NETWORK}/{PROTOCOL_DIRNAME_MAP[NETWORK][int(PROTOCOL)]}/specimens/ldms_feed/"
    fname = f"{NETWORK.lower()}.ldms{PROTOCOL}.{timestamp}.csv"
    # prev_fname = f"{NETWORK.lower()}.ldms{PROTOCOL}.{prev_timestamp}.csv"

    if os.path.exists(feed_dir + fname):
        new = pd.read_csv(feed_dir + fname, dtype=dtype_map)
        files = os.listdir(feed_dir)
        if len(files) > 1:
            prev_fname = np.sort(os.listdir(feed_dir))[-2]
            old = pd.read_csv(feed_dir + prev_fname, dtype=dtype_map)
        else:
            print(f"NO PRIOR LDMS SAVED FOR {NETWORK}{PROTOCOL}\n")
            old = pd.DataFrame()
    else:
        print(f"NO LDMS SAVED FOR {NETWORK}{PROTOCOL} TODAY\n")
        new = pd.DataFrame()

    # if os.path.exists(feed_dir + prev_fname):
    #     old = pd.read_csv(feed_dir + prev_fname, dtype=dtype_map)
    # else:
    #     print(f"NO LDMS SAVED FOR {NETWORK}{PROTOCOL} YESTERDAY\n")
    #     old = pd.DataFrame()
    # if os.path.exists(feed_dir + fname):
    #     new = pd.read_csv(feed_dir + fname, dtype=dtype_map)
    # else:
    #     print(f"NO LDMS SAVED FOR {NETWORK}{PROTOCOL} TODAY\n")
    #     new = pd.DataFrame()
    return old, new

def handle_affected_jobs(guspecs):
    handle_affected_jobs_old(guspecs)
    handle_affected_jobs_new(guspecs)

def handle_affected_jobs_old(guspecs): #report_outputs_corresponding_to_guspecs
    """
    given a list of guspecs with changes in ldms,
    report on any that are used by older (non-auto rerunnable) jobs
    """
    affected_outputs = set()
    for guspec in guspecs:
        try:
            for output in GUSPEC_TO_OUTPUT_PATH_OLD[guspec]:
                affected_outputs.add(output)
        except:
            print(f"{guspec} not in any results we're monitoring\n")
    if len(affected_outputs) > 0:
        print(f"AFFECTED OUTPUTS: {list(AFFECTED_OUTPUTS)}\n")
    else:
        print("NO PRE-2024 OUTPUTS AFFECTED\n")

def handle_affected_jobs_new(guspecs):
    """
    given a list guspecs with changes in ldms,
    return a list of yaml dicts corresponding to jobs that used that guspec
    """
    # pull all yaml_dicts associated with jobs
    print("finding yamls")
    yamls = [
        f for f in find_endpoints('/home/bhaddock/repos/sdmc-adhoc/processing_scripts', l=[]) if f.split("/")[-1]=="paths.yaml"
    ]
    yamls = [read_yaml(p) for p in yamls]
    affected_job_yamls = set()
    print("looping through guspecs, checking if in yamls")
    for guspec in guspecs:
        for y in yamls:
            if "guspecs" not in y.keys():
                print("adding guspecs to yaml")
                save_guspecs_to_yaml(yamls["yaml_path"])
                y = read_yaml(yamls["yaml_path"])
            if guspec in y["guspecs"]:
                print("found guspec in an output")
                affected_job_yamls.add(y)
    if len(affected_job_yamls) > 0:
        print("NEW OUTPUTS AFFECTED")
        for y in affected_job_yamls:
            output_path = get_output_path_from_yaml(y)
            print(f"The following output file affected: {output_path}")
            rerun_affected_job(y)
    else:
        print("NO NEW (2024+) OUTPUTS AFFECTED")

def rerun_affected_job(affected_job_yaml):
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

def get_output_path_from_yaml(yaml_dict):
    files = os.listdir(yaml_dict['savedir'])
    fname = np.sort([f for f in files if yaml_dict['output_prefix'] in f])[-1]
    output_path = yaml_dict['savedir']
    if output_path[-1] != "/":
        output_path += "/"
    output_path += fname
    return output_path

def find_endpoints(d, l):
    """
    for everything in dir ('/home/bhaddock/repos/sdmc-adhoc/processing_scripts')
    return if it's an endpoint (non-dir)
    """
    if os.path.isdir(d):
        for c in os.listdir(d):
            l = find_endpoints(d + "/" + c, l)
        return l
    else:
        return l + [d]

def read_yaml(yaml_path, add_path=True):
    """
    given yaml path
    optionally add yaml_path to yaml dict
    return yaml dict
    """
    with open(yaml_path, 'r') as file:
        yaml_dict = yaml.safe_load(file)
    yaml_dict["yaml_path"] = yaml_path
    return yaml_dict

def get_output_path_from_yaml(yaml_dict):
    files = os.listdir(yaml_dict['savedir'])
    fname = np.sort([f for f in files if yaml_dict['output_prefix'] in f])[-1]
    output_path = yaml_dict['savedir']
    if output_path[-1] != "/":
        output_path += "/"
    output_path += fname
    return output_path

def read_yaml(yaml_path, add_path=True):
    """
    given yaml path
    optionally add yaml_path to yaml dict
    return yaml dict
    """
    with open(yaml_path, 'r') as file:
        yaml_dict = yaml.safe_load(file)
    yaml_dict["yaml_path"] = yaml_path
    return yaml_dict

def save_guspecs_to_yaml(yaml_path):
    """
    given a filepath to a yaml that contains
    - a savedir
    - an output prefix(es)
    save the guspecs from the output to the yaml under key "guspecs"
    """
    yaml_dict = read_yaml(yaml_path, add_path=False)
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
            print(path)
            raise Exception(f"PUT EXCEPTION HERE; FILEPATH {path} WRONG")
        guspecs += data.guspec.tolist()
    yaml_dict['guspecs'] = list(set(guspecs))
    with open(yaml_path, 'w') as outfile:
        yaml.dump(yaml_dict, outfile, default_flow_style=False, sort_keys=False)
