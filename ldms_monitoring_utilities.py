import pandas as pd, numpy as np
import os
import sdmc_adhoc_processing.constants as constants
from ldms_monitoring_constants import *
import datetime
from typing import List



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
    N = 9 + len(str(PROTOCOL))
    feed_dir = f"/networks/vtn/lab/SDMC_labscience/studies/{NETWORK}/{PROTOCOL_DIRNAME_MAP[NETWORK][int(PROTOCOL)]}/specimens/ldms_feed/"
    things_in_dir = os.listdir(feed_dir)
    # sort from oldest to most recent, and only take the ones that look like hvtn.ldmsPROTOCOL
    applicable_things_in_dir = np.sort([i for i in things_in_dir if i[:N] == f'{NETWORK.lower()}.ldms{PROTOCOL}'])

    # delete all but the two most recent
    for file in applicable_things_in_dir[:-2]:
        print(f"DELETING THE FOLLOWING: {file}\n")
        os.remove(feed_dir + file)

def get_breakout_ldms_filepaths(PROTOCOL, NETWORK):
    today = datetime.date.today()
    timestamp = today.strftime("%Y%m%d")
    if today.weekday() == 0: # if it's a monday
        yesterday = datetime.date.today() - datetime.timedelta(days = 3)
    else: #any other day. note this is assuming the script isnt run sat or sun
        yesterday = datetime.date.today() - datetime.timedelta(days = 1)
    prev_timestamp = yesterday.strftime("%Y%m%d")

    feed_dir = f"/networks/vtn/lab/SDMC_labscience/studies/{NETWORK}/{PROTOCOL_DIRNAME_MAP[NETWORK][int(PROTOCOL)]}/specimens/ldms_feed/"
    fname = f"{NETWORK.lower()}.ldms{PROTOCOL}.{timestamp}.csv"
    prev_fname = f"{NETWORK.lower()}.ldms{PROTOCOL}.{prev_timestamp}.csv"

    if os.path.exists(feed_dir + prev_fname):
        old = pd.read_csv(feed_dir + prev_fname, dtype=dtype_map)
    else:
        print(f"NO LDMS SAVED FOR {NETWORK}{PROTOCOL} YESTERDAY\n")
        old = pd.DataFrame()
    if os.path.exists(feed_dir + fname):
        new = pd.read_csv(feed_dir + fname, dtype=dtype_map)
    else:
        print(f"NO LDMS SAVED FOR {NETWORK}{PROTOCOL} TODAY\n")
        new = pd.DataFrame()
    return old, new

def report_outputs_corresponding_to_guspecs(GUSPECS):
    AFFECTED_OUTPUTS = set()
    for guspec in GUSPECS:
        try:
            for output in GUSPEC_TO_OUTPUT_PATH[guspec]:
                AFFECTED_OUTPUTS.add(output)
        except:
            print(f"{guspec} not in any results we're monitoring\n")
    if len(AFFECTED_OUTPUTS) > 0:
        print(f"AFFECTED OUTPUTS: {list(AFFECTED_OUTPUTS)}\n")
    else:
        print("NO OUTPUTS AFFECTED\n")


def detect_ldms_diffs(PROTOCOL, NETWORK):
    """
    INPUT
        - protocol number (int)
        - network (HVTN or COVPN) associated with protocol
    FUNCTION: determines if yesterday / todays ldms has updates in any of the columns we're interested in
    """
    print(f"{NETWORK}{PROTOCOL}: reading in data")
    old, new = get_breakout_ldms_filepaths(PROTOCOL, NETWORK)

    # if either of these are empty, we have nothing to compare.
    if old.empty or new.empty:
        return

    print(f"{NETWORK}{PROTOCOL}: finding col diff")
    COL_DIFF = set(old.columns).symmetric_difference(new.columns)
    if len(COL_DIFF) > 0:
        print(f"COLUMNS DON'T MATCH. THE FOLLOWING NOT SHARED: {COL_DIFF}\n")

    print(f"{NETWORK}{PROTOCOL}: checking if lengths match")
    if len(old) != len(new):
        print(f"YESTERDAY AND TODAY DIFFERENT LENGTHS FOR {NETWORK} {PROTOCOL}; ATTEMPTING TO DE-DUP\n")
        old = old.drop_duplicates(ignore_index=True)
        new = new.drop_duplicates(ignore_index=True)
        print(f"{NETWORK}{PROTOCOL}: checking if lengths still don't match")
        if len(old) != len(new): # if lengths still dont match with deduping
            new_guspec_count = new.guspec.value_counts().to_frame()
            old_guspec_count = old.guspec.value_counts().to_frame()
            comp = new_guspec_count.merge(
                old_guspec_count,
                how='outer',
                left_index=True,
                right_index=True
            )
            DROPPED_GUSPECS = list(comp.loc[comp.count_x!=comp.count_y].index)
            print(f"{NETWORK}{PROTOCOL}: reporting on guspec diff:")
            if len(new) > len(old):
                print(f"THE FOLLOWING GUSPECS WERE ADDED TODAY: {DROPPED_GUSPECS}\n")
                report_outputs_corresponding_to_guspecs(DROPPED_GUSPECS)
                # TODO: LIST OUT STUDIES FROM THE APPLICABLE NETWORK/PROTOCOL
            if len(old) > len(new):
                print(f"THE FOLLOWING GUSPECS MISSING FROM NEW: {DROPPED_GUSPECS}\n")
                report_outputs_corresponding_to_guspecs(DROPPED_GUSPECS)
    # if columns and rows match, ensure sorting is the same
    print(f"{NETWORK}{PROTOCOL}: checking if lengths now match:")
    if len(old) == len(new):
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

            print(f"{NETWORK}{PROTOCOL}: storing dataframe diff")
            # catch any rows that have mismatch
            diff = new.compare(old)
            if len(diff)>0:
                print(f"CHANGES DETECTED IN EXISTING GUSPECS\n")
                AFFECTED_GUSPECS = list(new.iloc[diff.index].guspec.unique())
                print(f"AFFECTED GUSPECS: {AFFECTED_GUSPECS}\n")
                report_outputs_corresponding_to_guspecs(AFFECTED_GUSPECS)
            else:
                print(f"NO DIFF DETECTED FOR {NETWORK}{PROTOCOL}\n")
