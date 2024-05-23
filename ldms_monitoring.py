from ldms_monitoring_utilities import *

yamlpath = os.path.dirname(__file__) + "/constants.yaml"
with open(yamlpath, 'r') as file:
    yamldict = yaml.safe_load(file)

hvtn_protocols = yamldict["hvtn_protocols"]
covpn_protocols = yamldict["covpn_protocols"]

def check_for_updates(protocols: list, network: str) -> None:
    # get today's timestamp
    updated = was_ldms_updated(network)
    if updated:
        # save new copies of ldms
        if len(protocols)>0:
            # if there are at least two ldmss saved, delete the older ones
            save_todays_ldms(protocols, network)
        #     for protocol in protocols:
        #         delete_old_ldms(protocol, network)

        # determine if there are any updates, output corresponding alerts
        for protocol in protocols:
            try:
                detect_ldms_diffs(int(protocol), network)
            except:
                print(f"ISSUE WITH {network}{int(protocol)}; CHECK CODE")
    else:
        print(f"{network} LDMS HASN'T BEEN UPDATED TODAY\n")

check_for_updates(hvtn_protocols, "HVTN")
check_for_updates(covpn_protocols, "CoVPN")
