import pandas as pd

# paths to outputs we're monitoring ------------------------------------------##
OUTPUT_REGISTRY_OLD = [
 '/networks/vtn/lab/SDMC_labscience/studies/CoVPN/CoVPN3008/assays/ics/misc_files/CoVPN3008_FCM08_AP59_CH_20221229E01FLO_B.txt',
 '/networks/vtn/lab/SDMC_labscience/studies/CoVPN/CoVPN3008/assays/ics/misc_files/CoVPN3008_FCM08_AP59_CH_20221229E01FLO_C.txt',
 '/networks/vtn/lab/SDMC_labscience/studies/CoVPN/CoVPN3008/assays/ics/misc_files/summary_ds_2023-04-24.csv',
 '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN Auxiliary Studies/111-ANS_Maenetje_HVTN_503/draft_FCM08_output_20210616.txt',
 '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN Auxiliary Studies/Mycoplasma/161-EXS_Mycoplasma_Mayo_HVTN_092_096/sdmc_output/161EXS_Mycoplasma_Mayo_HVTN_092_096_20190312.txt',
 '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN Auxiliary Studies/Mycoplasma/161-EXS_Mycoplasma_Mayo_HVTN_092_096/sdmc_output/161EXS_Mycoplasma_Mayo_HVTN_092_096_20190313.txt',
 '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN Auxiliary Studies/Mycoplasma/162-EXS_Mycoplasma_Iowa_HVTN_092_096/sdmc_output/162EXS_Mycoplasma_Iowa_HVTN_092_096_Part1_20190225.txt',
 '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN Auxiliary Studies/Mycoplasma/162-EXS_Mycoplasma_Iowa_HVTN_092_096/sdmc_output/162EXS_Mycoplasma_Iowa_HVTN_092_096_Part3_20190711.txt',
 '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN Auxiliary Studies/Mycoplasma/162-EXS_Mycoplasma_Iowa_HVTN_092_096/sdmc_output/162EXS_Mycoplasma_Iowa_HVTN_092_096_Part3_20191014.txt',
 '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN139/assays/NAb/Wistar/misc_files/HVTN139_Adenovirus_nAb_Wistar_2023-08-30.csv',
 '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN139/assays/NAb/Wistar/misc_files/archive/HVTN139_Adenovirus_nAb_Wistar_2023-07-17.csv',
 '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN302/assays/AE_assays/basophil_activation_testing/misc_files/HVTN302_basophil_activation_processed_2023-11-22.txt',
 '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN303/assays/AE_assays/Circulating Immune Complex Panel/misc_files/data_processing/DRAFT_HVTN303_ARUP_C1C-C1q_processed_2023-09-25.txt',
 '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN303/assays/MSD/Binding_antibody_via_MSD/misc_files/data_processing/DRAFT_HVTN303_MSD_binding_processed_2023-09-26.txt',
 '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN303/assays/MSD/Cytokines_via_MSD/misc_files/data_processing/HVTN303_MSD_cytokines_processed_2024-02-16.txt',
]

def find_endpoints(d, l):
    if os.path.isdir(d):
        for c in os.listdir(d):
            l = find_endpoints(d + "/" + c, l)
        return l
    else:
        return l + [d]

def get_output_filepath_from_yaml(yamlpath: str) -> list:
    with open(yamlpath, 'r') as file:
        yaml_dict = yaml.safe_load(file)
    output_dir = yaml_dict["savedir"]
    if output_dir[-1]!="/":
        output_dir += "/"
    #TODO: add try/excepts for if there's a nonmatch or this is missing
    try:
        prefixes = yaml_dict["output_prefix"]
        if not isinstance(prefixes, list):
            prefixes = [prefixes]
        fnames = []
        for prefix in prefixes:
            matches = [file for file in os.listdir(output_dir) if file[:len(prefix)]==prefix]
            fnames += [output_dir + np.sort(matches)[-1]]
        return fnames
    except:
        return [output_dir + "MISSING"]

def get_output_registry_new():
    endpoints = find_endpoints('/home/bhaddock/repos/sdmc-adhoc/processing_scripts', l=[])
    yamls = [file for file in endpoints if file.split("/")[-1]=="paths.yaml"]
    output_paths = []
    for y in yamls:
        output_paths += get_output_filepath_from_yaml(y)
    return output_paths

OUTPUT_REGISTRY = get_output_registry_new()

# list of outputs we're monitoring -------------------------------------------##
def read_in_output(path):
    """
    INPUT: filepath to csv or txt of data
    OUTPUT: data contained at filepath, columns standardized
    """
    filetype = path.rpartition(".")[-1]
    if filetype =="txt":
        data = pd.read_csv(path, sep="\t")
    elif filetype == "csv":
        data = pd.read_csv(path)
    else:
        return f"UNRECOGNIZED FILETYPE: {filetype}"
    # if we were able to read in the data
    data.columns = [i.lower() for i in data.columns]
    return data

OUTPUTS_OLD = [read_in_output(path) for path in OUTPUT_REGISTRY_OLD]
OUTPUTS = [read_in_output(path) for path in OUTPUT_REGISTRY]

PATH_OUTPUT_PAIRS = zip(OUTPUT_REGISTRY + OUTPUT_REGISTRY_OLD, ALL_OUTPUTS + ALL_OUTPUTS_OLD)

# DTYPES FOR LDMS ------------------------------------------------------------##
dtype_map = {
    "txtpid":str,
    "drawdm":'Int64',
    "drawdd":'Int64',
    "drawdy":'Int64',
    "vidval":'Float64',
    "lstudy":'Float64',
    "guspec":str,
    "primstr":str,
    "addstr":str,
    "dervstr":str
    }

# list of corresponding guspecs ----------------------------------------------##
def get_protocol_gus(output):
    """
    For a given output, returns all guspecs associated with the output's
    protocol(s)
    """
    if "guspec" in output.columns:
        d = {}
        for protocol in output.protocol.unique():
            d[protocol] = list(output.loc[output.protocol==protocol].guspec.unique())
        return d
    else:
        return ["NO GUSPEC COLUMN"]

GUSPEC_PROTOCOL_PER_OUTPUT = {
    path:get_protocol_gus(output) for (path,output) in PATH_OUTPUT_PAIRS
}

# list of outputs corresponding to guspecs -----------------------------------##
def generate_guspec_to_path_map(PATH_OUTPUT_PAIRS):
    bigmap = {}
    for (path,output) in PATH_OUTPUT_PAIRS:
        if "guspec" in output.columns:
            for guspec in output.guspec.unique():
                if guspec not in bigmap:
                    bigmap[guspec] = [path]
                else:
                    bigmap[guspec] += [path]

    return bigmap

GUSPEC_TO_OUTPUT_PATH_OLD = generate_guspec_to_path_map(zip(OUTPUT_REGISTRY_OLD, ALL_OUTPUTS_OLD))

# dict of protocol paths for saving ldms -------------------------------------##
PROTOCOL_DIRNAME_MAP = {"HVTN": {49: 'HVTN049',
                                84: 'HVTN084',
                                87: 'HVTN087',
                                88: 'HVTN088',
                                90: 'HVTN090',
                                94: 'HVTN094',
                                97: 'HVTN097',
                                98: 'HVTN098',
                                100: 'HVTN100',
                                105: 'HVTN105',
                                106: 'HVTN106',
                                107: 'HVTN107',
                                108: 'HVTN108',
                                111: 'HVTN111',
                                114: 'HVTN114',
                                115: 'HVTN115',
                                116: 'HVTN116',
                                117: 'HVTN117',
                                118: 'HVTN118',
                                119: 'HVTN119',
                                120: 'HVTN120',
                                121: 'HVTN121',
                                122: 'HVTN122',
                                123: 'HVTN123',
                                124: 'HVTN124',
                                125: 'HVTN125',
                                127: 'HVTN127_HPTN087',
                                128: 'HVTN128',
                                129: 'HVTN129_HPTN088',
                                130: 'HVTN130_HPTN089',
                                131: 'HVTN131_HPTN090',
                                132: 'HVTN132',
                                133: 'HVTN133',
                                134: 'HVTN134',
                                135: 'HVTN135',
                                136: 'HVTN136_HPTN092',
                                137: 'HVTN137',
                                139: 'HVTN139',
                                140: 'HVTN140_HTPN101',
                                142: 'HVTN142',
                                144: 'HVTN144',
                                300: 'HVTN300',
                                301: 'HVTN301',
                                302: 'HVTN302',
                                303: 'HVTN303',
                                304: 'HVTN304',
                                305: 'HVTN305',
                                306: 'HVTN306',
                                307: 'HVTN307',
                                308: 'HVTN308',
                                310: 'HVTN310',
                                311: 'HVTN311',
                                405: 'HVTN405_HPTN1901',
                                502: 'HVTN502',
                                503: 'HVTN503',
                                505: 'HVTN505',
                                602: 'HVTN602',
                                603: 'HVTN603_ACTG5397',
                                604: 'HVTN604',
                                605: 'HVTN605',
                                702: 'HVTN702',
                                703: 'HVTN703_704',
                                704: 'HVTN704',
                                705: 'HVTN705',
                                706: 'HVTN706',
                                804: 'HVTN804_HPTN095',
                                805: 'HVTN805_HPTN093',
                                806: 'HVTN806_HPTN108_A5416',
                                807: 'HVTN807',
                                92: 'HVTN92',
                                96: 'HVTN96'},
                                "CoVPN": {3001: 'CoVPN3001',
                                         3002: 'CoVPN3002',
                                         3003: 'CoVPN3003',
                                         3004: 'CoVPN3004',
                                         3005: 'CoVPN3005',
                                         3006: 'CoVPN3006',
                                         3007: 'CoVPN3007',
                                         3008: 'CoVPN3008',
                                         5001: 'CoVPN5001',
                                         5002: 'CoVPN5002',
                                         5003: 'CoVPN5003'}
                                }
