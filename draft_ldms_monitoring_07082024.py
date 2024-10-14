# ---------------------------------------------------- #
# LDMS MONITORING
# 07-08-2024
# Author: Beatrix Haddock
# Every day, for specified protocols:
#   - save copy of ldms
#   - delete copies of ldms older than X
#   - check if ldms has changed between today
#   prev saved version.
#   - send email notification of any relevant changes
# -----------------------------------------------------#



def monitor_ldms(protocol_list: list, network: str) -> None:
    pass


if __name__=="__main__":
    # set up reporting
    CHANGES_TO_REPORT = ""
    PROTOCOLS_WITHOUT_CHANGE = ""
    
    # get list of protocols for hvtn and covpn that we're tracking
    yamlpath = os.path.dirname(__file__) + "/constants.yaml"
    with open(yamlpath, 'r') as file:
        yamldict = yaml.safe_load(file)

    hvtn_protocols = yamldict["hvtn_protocols"]
    covpn_protocols = yamldict["covpn_protocols"]

    # report which protocols we're monitoring
    print("Monitoring the following protocols:")
    print(f"HVTN: {", ".join([str(p) for p in hvtn_protocols])}\n")
    print(f"CoVPN: {", ".join([str(p) for p in covpn_protocols])}\n")

    # check ldms for changes
    monitor_ldms(hvtn_protocols, "HVTN")
    monitor_ldms(covpn_protocols, "CoVPN")
