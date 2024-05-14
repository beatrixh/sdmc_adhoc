from ldms_monitoring_utilities import *
from ldms_monitoring_constants import *

# dataframe of network-protocol pairs in outputs we're monitoring
NETWORK_PROTOCOLS = pd.concat([data[["network", "protocol"]].drop_duplicates() for data in ALL_OUTPUTS]).dropna()
NETWORK_PROTOCOLS.protocol = NETWORK_PROTOCOLS.protocol.astype(int).astype(float)
NETWORK_PROTOCOLS = NETWORK_PROTOCOLS.drop_duplicates()

# protocols associated with hvtn vs covpn
HVTN_PROTOCOLS = list(NETWORK_PROTOCOLS.loc[NETWORK_PROTOCOLS.network=="HVTN"].protocol)
COVPN_PROTOCOLS = list(NETWORK_PROTOCOLS.loc[NETWORK_PROTOCOLS.network=="CoVPN"].protocol)

def check_for_updates(PROTOCOLS, NETWORK):
    # get today's timestamp
    updated = was_ldms_updated(NETWORK)
    if updated:
        # save new copies of ldms
        if len(PROTOCOLS)>0:
            # if there are at least two ldmss saved, delete the older ones
            save_todays_ldms(PROTOCOLS, NETWORK)
            for protocol in PROTOCOLS:
                delete_old_ldms(protocol, NETWORK)

        # determine if there are any updates, output corresponding alerts
        for protocol in PROTOCOLS:
            try:
                detect_ldms_diffs(int(protocol), NETWORK)
            except:
                print(f"ISSUE WITH {NETWORK}{int(protocol)}; CHECK CODE")
    else:
        print(f"{NETWORK} LDMS HASN'T BEEN UPDATED TODAY\n")

check_for_updates(HVTN_PROTOCOLS, "HVTN")
check_for_updates(COVPN_PROTOCOLS, "CoVPN")

## ---------------------------------------------------------------------------##
#TODO: store dictionaries so don't have to read inputs in every time
# primaryMap = {path -> network -> protocol -> guspec}
# networkProtocolMap = {network: protocols}
# trackingMap = {path -> last updated}

# def stayUpToDate():
#     if current_last_updated > last_updated:
#         update primaryMap
#         update trackingMap
#         update networkProtocolMap
#
# def checkForUpdates():
#
#
# if ldms updated today:
#     - run stayUpToDate
#     - saveTodaysLdms(protocols = networkProtocolmap[network], network = network)
#     - delete older ldmss
#     - check for updates
# else:
#     notify in logs that rerun needed
#     #todo: rerun automatically
## ---------------------------------------------------------------------------##
