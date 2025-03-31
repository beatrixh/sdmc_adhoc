## ---------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 03/12/2025
# Purpose:
#     - Save back-ups of subsets of LDMS
#     - Check for changes since yesterday
## ---------------------------------------------------------------------------##
import pandas as pd, numpy as np
import os
import yaml
import datetime
from typing import List
import sdmc_tools.constants as constants
from yaml_handling import *
from access_ldms import *

def main():
    protocols = {}
    protocols['hvtn'] = yamldict["hvtn_protocols"]
    protocols['covpn'] = yamldict["covpn_protocols"]

    for network in ['hvtn','covpn']:
        for protocol in protocols[network]:
            save_todays_ldms(protocol, network)
            delete_old_ldms(protocol, network)
            detect_ldms_diffs(protocol, network)

## constants and functions -------------------------------------------------- ##
yamlpath = os.path.dirname(__file__) + "/constants.yaml"
with open(yamlpath, 'r') as file:
    yamldict = yaml.safe_load(file)

protocols = {}
protocols['hvtn'] = yamldict["hvtn_protocols"]
protocols['covpn'] = yamldict["covpn_protocols"]

PROTOCOL_DIRNAME_MAP = yamldict["PROTOCOL_DIRNAME_MAP"]
GUSPEC_TO_OUTPUT_PATH_OLD = yamldict["GUSPEC_TO_OUTPUT_PATH_OLD"]

def detect_ldms_diffs(PROTOCOL, NETWORK):
    if NETWORK.lower() == 'covpn':
        NETWORK = 'CoVPN'
    elif NETWORK.lower() == 'hvtn':
        NETWORK = 'HVTN'
    else:
        raise Exception(f"NETWORK must be one of 'HVTN' or 'CoVPN' (case insensitive). Submitted {NETWORK}")

    ## pull old and new LDMS ------------------------------------------------ ##
    today = datetime.date.today()
    timestamp = today.strftime("%Y%m%d")

    feed_dir = f"/networks/vtn/lab/SDMC_labscience/studies/{NETWORK}/{PROTOCOL_DIRNAME_MAP[NETWORK.upper()][int(PROTOCOL)]}/specimens/ldms_feed/"
    fname = f"{NETWORK.lower()}.ldms{PROTOCOL}.{timestamp}.csv"

    ## save LDMS for today if it doesn't exist
    if not os.path.exists(feed_dir + fname):
        save_todays_ldms(PROTOCOL, NETWORK)

    ## read in todays LDMS
    new = pd.read_csv(feed_dir + fname, dtype=constants.LDMS_DTYPE_MAP)
    ## read in previous LDMS
    files = os.listdir(feed_dir)
    if len(files) > 1:
        files = os.listdir(feed_dir)
        files = [i for i in files if f'{NETWORK.lower()}.ldms{PROTOCOL}' in i]
        prev_fname = np.sort(files)[-2]
        old = pd.read_csv(feed_dir + prev_fname, dtype=constants.LDMS_DTYPE_MAP)
        # print(f"{NETWORK}{PROTOCOL}: Comparing {fname} and {prev_fname}")
    else:
        print(f"NO PRIOR LDMS SAVED FOR {NETWORK}{PROTOCOL}")
        old = pd.DataFrame()

    # if either of these are empty, we have nothing to compare ----------------#
    if old.empty or new.empty:
        return

    # the columns have to match -- otherwise have big problems ----------------#
    # print(f"{NETWORK}{PROTOCOL}: Checking for col diff")
    COL_DIFF = set(old.columns).symmetric_difference(new.columns)
    if len(COL_DIFF) > 0:
        print(f"{NETWORK}{PROTOCOL}: COLUMNS DON'T MATCH. THE FOLLOWING NOT SHARED: {COL_DIFF}")
        return

    # ensure column types / sorting match -------------------------------------#
    col_sort_order = ['guspec', 'txtpid', 'drawdm', 'drawdd', 'drawdy', 'vidval', 'lstudy',
       'primstr', 'addstr', 'dervstr']

    for df in [old, new]:
        for col in ['drawdm','drawdd','drawdy']:
            df[col] = df[col].astype(int)
        for col in ['vidval','lstudy']:
            df[col] = df[col].astype(float)
        df['txtpid'] = df.txtpid.astype(str)
        df = df.sort_values(by=col_sort_order, ignore_index=True).drop_duplicates()

    # were any guspecs removed ------------------------------------------------#
    guspecs_removed = set(old.guspec).difference(new.guspec)
    if guspecs_removed:
        print(f"{NETWORK}{PROTOCOL}: THE FOLLOWING GUSPECS HAVE BEEN REMOVED FROM LDMS: {guspecs_removed}")
        handle_affected_jobs(guspecs_removed)

    # report on any added
    guspecs_added = set(new.guspec).difference(old.guspec)
    if guspecs_added:
        print(f"{NETWORK}{PROTOCOL}: THE FOLLOWING GUSPECS WERE ADDED TO LDMS: {guspecs_added}")

    # of the guspecs in both, is there still the name number of rows per guspec
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

    # did anything change among the guspecs that are in both ------------------#
    # sort for comparison
    new = new.sort_values(by=col_sort_order, ignore_index=True)
    old = old.sort_values(by=col_sort_order, ignore_index=True)

    # compare
    diff = new.replace("N/A", np.nan).compare(old.replace("N/A", np.nan))
    if len(diff) > 0:
        print(f"CHANGES DETECTED IN SHARED GUSPECS")
        AFFECTED_GUSPECS = list(new.iloc[diff.index].guspec.unique())
        print(f"AFFECTED GUSPECS: {AFFECTED_GUSPECS}")
        handle_affected_jobs(AFFECTED_GUSPECS)

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
    # print("Checking to see if new (2024+) outputs affected\n")
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
            print(f"! The following output file affected: {output_path}")
            # rerun_affected_job(y)
    else:
        print("\nNO NEW (2024+) OUTPUTS AFFECTED")

## delete old ldms ---------------------------------------------------------- ##
def delete_old_ldms(PROTOCOL: str, NETWORK: str) -> None:
    if NETWORK.lower() == 'covpn':
        NETWORK = 'CoVPN'
    elif NETWORK.lower() == 'hvtn':
        NETWORK = 'HVTN'
    else:
        raise Exception(f"NETWORK must be one of 'HVTN' or 'CoVPN' (case insensitive). Submitted {NETWORK}")
    PROTOCOL = int(PROTOCOL)
    N = len(NETWORK) + len(".ldms") + len(str(PROTOCOL))
    feed_dir = f"/networks/vtn/lab/SDMC_labscience/studies/{NETWORK}/{PROTOCOL_DIRNAME_MAP[NETWORK.upper()][int(PROTOCOL)]}/specimens/ldms_feed/"
    things_in_dir = os.listdir(feed_dir)
    # sort from oldest to most recent, and only take the ones that look like hvtn.ldmsPROTOCOL
    applicable_things_in_dir = np.sort([i for i in things_in_dir if i[:N] == f'{NETWORK.lower()}.ldms{PROTOCOL}'])[::-1]

    # delete all but the thirty most recent
    for file in applicable_things_in_dir[30:]:
        # print(f"DELETING THE FOLLOWING: {file}\n")
        os.remove(feed_dir + file)

def save_todays_ldms(PROTOCOL: str, NETWORK: str) -> None:
    """
    INPUTS
     - PROTOCOLS: list of PROTOCOLS (numeric) from one network
     - NETWORK: string ("HVTN" or "CoVPN") corresponding to above PROTOCOLS
    """
    if NETWORK.lower() == 'covpn':
        NETWORK = 'CoVPN'
    elif NETWORK.lower() == 'hvtn':
        NETWORK = 'HVTN'
    else:
        raise Exception(f"NETWORK must be one of 'HVTN' or 'CoVPN' (case insensitive). Submitted {NETWORK}")

    # read in ldms
    ldms = pull_one_protocol(NETWORK, PROTOCOL)

    # save ldms
    timestamp = datetime.date.today().strftime("%Y%m%d")
    network_path = f"/networks/vtn/lab/SDMC_labscience/studies/{NETWORK}/"
    protocol_path = network_path + PROTOCOL_DIRNAME_MAP[NETWORK.upper()][int(PROTOCOL)]
    if not os.path.exists(protocol_path):
        print(f"creating protocol dir for {int(PROTOCOL)}")
        os.mkdir(protocol_path)
    specimen_path = protocol_path + "/specimens"
    if not os.path.exists(specimen_path):
        print(f"creating specimens dir for {int(PROTOCOL)}")
        os.mkdir(specimen_path)
    ldms_feed_path = specimen_path + "/ldms_feed"
    if not os.path.exists(ldms_feed_path):
        print(f"creating ldms_feed dir for {int(PROTOCOL)}")
        print(ldms_feed_path)
        os.mkdir(ldms_feed_path)

    savepath = ldms_feed_path + f"/{NETWORK.lower()}.ldms{int(PROTOCOL)}.{timestamp}.csv"
    if not os.path.exists(savepath):
        # print(f"CREATING THE FOLLOWING: {savepath}")
        ldms.to_csv(savepath, index=False)
    else:
        print(f"{savepath} already exists")


if __name__=="__main__":
    main()
