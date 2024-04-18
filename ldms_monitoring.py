from ldms_monitoring_utilities import *
from ldms_monitoring_constants import *

# dataframe of network-protocol pairs in outputs we're monitoring
NETWORK_PROTOCOLS = pd.concat([data[["network", "protocol"]].drop_duplicates() for data in ALL_OUTPUTS]).dropna()
NETWORK_PROTOCOLS.protocol = NETWORK_PROTOCOLS.protocol.astype(int).astype(float)
NETWORK_PROTOCOLS = NETWORK_PROTOCOLS.drop_duplicates()

# protocols associated with hvtn vs covpn
HVTN_PROTOCOLS = list(NETWORK_PROTOCOLS.loc[NETWORK_PROTOCOLS.network=="HVTN"].protocol)
COVPN_PROTOCOLS = list(NETWORK_PROTOCOLS.loc[NETWORK_PROTOCOLS.network=="CoVPN"].protocol)

# get timestamp that LDMS most recently updated
LAST_TOUCHED = os.stat(constants.LDMS_PATH_HVTN).st_mtime
LAST_TOUCHED = datetime.datetime.fromtimestamp(LAST_TOUCHED)
#TODO: condition on COVPN LAST TOUCHED
# if LDMS has been updated today
if LAST_TOUCHED.date() == datetime.date.today():
    # save new copies of ldms
    if len(HVTN_PROTOCOLS)>0:
        save_todays_ldms(HVTN_PROTOCOLS, "HVTN")
        # if there are at least two ldmss saved, delete the older ones
        # for PROTOCOL in HVTN_PROTOCOLS:
        #     delete_old_ldms(PROTOCOL, "HVTN")
    if len(COVPN_PROTOCOLS)>0:
        save_todays_ldms(COVPN_PROTOCOLS, "CoVPN")
        # if there are at least two ldmss saved, delete the older ones
        # for PROTOCOL in COVPN_PROTOCOLS:
        #     delete_old_ldms(PROTOCOL, "CoVPN")
    # determine if there are any updates, output corresponding alerts
    for protocol in HVTN_PROTOCOLS:
        detect_ldms_diffs(int(protocol), "HVTN")
    for protocol in COVPN_PROTOCOLS:
        detect_ldms_diffs(int(protocol), "CoVPN")
elif LAST_TOUCHED.date() < datetime.date.today():
    print("LDMS HASN'T BEEN UPDATED TODAY\n")


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
