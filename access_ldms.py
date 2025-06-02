import pandas as pd
import numpy as np
import psycopg
import yaml

STANDARD_COLS = ['txtpid',
 'drawdm',
 'drawdd',
 'drawdy',
 'vidval',
 'lstudy',
 'guspec',
 'primstr',
 'addstr',
 'dervstr']

# grab my password
yaml_path = "/home/bhaddock/repos/config.yaml"
with open(yaml_path, 'r') as file:
    config = yaml.safe_load(file)

def get_cursor():
    """
    Get cursor into datamart
    """
    conn = psycopg.connect(host="datamart.scharp.org",
                            dbname="datamart02_realtime",
                            user="bhaddock",
                            password=config['password'])

    return conn.cursor()

def pull_protocol_map() -> dict:
    """
    Return a dictionary that maps each protocol to the schema it is stored in.
    e.g.,
    {'hvtn82': 'hvtn082',
    'hvtn139': 'hvtn139_dummy_arm',
    'hvtn503': 'hvtn503',
    'covpn3002': 'covpn3002'}
    -----
    DICTIONARY ARGS
    - keys are of the form {network}{protocol}
        - with no spaces,
        - where 'network' is all lowercase,
        - protocol is int-formatted with no leading zeros.
        - EXCEPT: 'unmapped' and 'hvtn_unmapped' map to the unmapped hvtn protocol schema
    """
    # pull all schemas that have an imported_fstrf_specimen_aliquot table---------------------------------------##
    cursor = get_cursor()
    cursor.execute("""SELECT table_schema,table_name FROM information_schema.tables where table_name like 'imported_fstrf_specimen_aliquot'""")
    tables = cursor.fetchall()

    schemas = np.array(tables)[:,0]

    schema_usemap = {}
    for s in schemas:
        protocol = s.split("_")[0]
        if protocol not in schema_usemap:
            schema_usemap[protocol] = [s]
        else:
            schema_usemap[protocol] += [s]

    # grab all protocols with only one schema per protocol -----------------------------------------------------##
    one_option_only = pd.DataFrame({k:v for (k,v) in schema_usemap.items() if len(v) == 1}).T.reset_index()
    one_option_only.columns = ['network_protocol','Schema']
    one_option_only['Use'] = True

    one_option_only.loc[one_option_only.Schema=='hvtn_unmapped','network_protocol'] = 'hvtn_unmapped'

    # now deal with protocols that have multiple schemas
    multiples = {k:v for (k,v) in schema_usemap.items() if len(v) > 1}

    which_to_use = pd.DataFrame([i for j in multiples.values() for i in j], columns=['schema'])
    which_to_use['use'] = False
    which_to_use = pd.concat([which_to_use.schema.str.split("_", expand=True)[0], which_to_use], axis=1)
    which_to_use.columns = ['network_protocol', 'Schema', 'Use']

    which_to_use.loc[which_to_use.Schema=="hvtn144_dummy_arm", 'Use'] = True
    which_to_use.loc[which_to_use.Schema=="hvtn503", 'Use'] = True
    which_to_use.loc[which_to_use.Schema=="hvtn139_dummy_arm", 'Use'] = True
    which_to_use.loc[which_to_use.Schema=="hvtn312_dummy_arm", 'Use'] = True
    which_to_use.loc[which_to_use.Schema=="hvtn605_dummy_arm", 'Use'] = True
    which_to_use.loc[which_to_use.Schema=="hvtn320_dummy_arm", 'Use'] = True

    # # add in HVTN503.1
    # which_to_use.loc[len(which_to_use)] = ['hvtn503.1', 'hvtn503_1', True]
    # which_to_use.loc[len(which_to_use)] = ['hvtn503_1', 'hvtn503_1', True]

    # ensure that for all protocols with multiple schemas, we chose exactly one
    check_one_per_protocol = which_to_use.groupby('network_protocol').Use.sum()
    issues = check_one_per_protocol.loc[check_one_per_protocol != 1].index.tolist()
    if len(issues) > 0:
        raise Warning(f"For the following protocols, there are multiple schemas! Need to choose which to use: {issues}")

    which_to_use = pd.concat([which_to_use, one_option_only], axis=0)

    ## remove leading zeros in network_protocol names ---------------------------------------------------------##

    # drop the unmapped network so can cast others to int
    which_to_use = which_to_use.loc[~which_to_use.network_protocol.isin(['hvtn_unmapped', 'hvtn503_1', 'hvtn503.1'])]

    # pull out network
    which_to_use['network'] = ['hvtn' if 'hvtn' in i else 'covpn' for i in which_to_use.network_protocol]

    # pull out network
    which_to_use['protocol'] = [int(i[4:]) if i[:4]=='hvtn' else int(i[5:]) for i in which_to_use.network_protocol]

    # concat network/protocol back together
    which_to_use['key'] = which_to_use.network + which_to_use.protocol.astype(str)

    # convert to dictionary, add unmapped back in
    # WILL NEED TO UPDATE IF THEY ADD A COVPN UNMAPPED
    which_to_use = which_to_use.loc[which_to_use.Use].set_index('key')['Schema'].to_dict()
    which_to_use['unmapped'] = 'hvtn_unmapped'
    which_to_use['hvtn_unmapped'] = 'hvtn_unmapped'
    which_to_use['hvtn503.1'] = 'hvtn503_01'
    which_to_use['hvtn503_1'] = 'hvtn503_01'

    return which_to_use

def pull_one_protocol(network, protocol, usecols=STANDARD_COLS):
    """
    ----
    INPUT
    network: network name. one of 'hvtn' or 'covpn'

    usecols
    -----
    TODO: look into proper way to call conn.close()
    """
    if protocol=='unmapped':
        network_protocol = 'hvtn_unmapped'
    elif protocol in ['503.1', '503_1', 503.1, '503_01', '503.01', 503.01]:
        network_protocol = 'hvtn503_01'
    else:
        if isinstance(protocol, str):
            protocol = protocol.lower().replace(" ","")
        network = network.lower().replace(" ","").replace("_","")
        pmap = pull_protocol_map()
        schemas = pd.DataFrame({'network_protocol':pmap.keys(), 'schema_name':pmap.values()})
        schemas['network'] = ['hvtn' if 'hvtn' in i else 'covpn' for i in schemas.network_protocol]
        schemas['protocol'] = [i[4:] if i[:4]=='hvtn' else i[5:] for i in schemas.network_protocol]
        schemas.loc[schemas.schema_name.str.contains("unmapped"), 'protocol'] = np.nan
        schemas['protocol'] = schemas.protocol.astype(float)

        schema = schemas.loc[(
            (schemas.network==network) &
            (schemas.protocol==float(protocol))
        )]
        if len(schema) > 1:
            raise Warning("ERROR: Network/Protocol combo returned more than one table")

        network_protocol = f"{network}{int(protocol)}"
        network_protocol = pmap[network_protocol]

    # stringify columns to pull
    if isinstance(usecols, list):
        columns_to_pull = ', '.join(usecols)
    # pull all if user passes 'all'
    elif usecols == 'all':
        columns_to_pull = "*"
    else:
        columns_to_pull = usecols

    # get db cursor
    cursor = get_cursor()

    # pull data
    cursor.execute(f"""SELECT {columns_to_pull} FROM {network_protocol}.imported_fstrf_specimen_aliquot""")
    data = cursor.fetchall()

    # if pulled all columns, get them as a list
    if usecols == 'all':
        cursor.execute(f"""select column_name from information_schema.columns where table_schema='{network_protocol}' AND table_name like 'imported_fstrf_specimen_aliquot' order by ordinal_position;""")
        usecols = cursor.fetchall()
        usecols = [i[0] for i in usecols]
    if not isinstance(usecols, list):
        usecols = [usecols]

    # close connection
    cursor.close()

    return pd.DataFrame(data, columns=usecols)

def pull_multiple_protocols(network, protocols, usecols=STANDARD_COLS):

    # stringify columns to pull
    if isinstance(usecols, list):
        columns_to_pull = ', '.join(usecols)
    # pull all if user passes 'all'
    elif usecols == 'all':
        columns_to_pull = "*"
    else:
        columns_to_pull = usecols

    # format network_protocols correctly
    network = network.lower().replace(" ","").replace("_","")
    pmap = pull_protocol_map()
    network_protocols = []
    if protocols=='all':
        network_protocols = list(set(pmap.values()))
        network_protocols = [i for i in network_protocols if network in i]
    else:
        for protocol in protocols:
            if protocol=='unmapped':
                network_protocols += ['hvtn_unmapped']
            else:
                if isinstance(protocol, str):
                    protocol = protocol.lower().replace(" ","")
                schemas = pd.DataFrame({'network_protocol':pmap.keys(), 'schema_name':pmap.values()})
                schemas['network'] = ['hvtn' if 'hvtn' in i else 'covpn' for i in schemas.network_protocol]
                schemas['protocol'] = [i[4:] if i[:4]=='hvtn' else i[5:] for i in schemas.network_protocol]
                schemas.loc[schemas.schema_name.str.contains("unmapped"), 'protocol'] = np.nan
                schemas['protocol'] = schemas.protocol.astype(float)

                schema = schemas.loc[(
                    (schemas.network==network) &
                    (schemas.protocol==float(protocol))
                )]
                if len(schema) > 1:
                    raise Warning(f"ERROR: {network}/{protocol} combo returned more than one table")

                network_protocol = f"{network}{int(protocol)}"
                network_protocol = pmap[network_protocol]
                network_protocols += [network_protocol]

    # connect to db
    cursor = get_cursor()

    # build query
    query_string = ''
    for p in network_protocols:
        query_string += f' UNION ALL SELECT {columns_to_pull} FROM {p}.imported_fstrf_specimen_aliquot'

    query_string = query_string[11:]

    # pull data
    cursor.execute(query_string)
    data = cursor.fetchall()

    # if pulled all columns, get them as a list
    if usecols == 'all':
        cursor.execute(f"""select column_name from information_schema.columns where table_schema='{network_protocols[0]}' AND table_name like 'imported_fstrf_specimen_aliquot' order by ordinal_position;""")
        usecols = cursor.fetchall()
        usecols = [i[0] for i in usecols]
    if not isinstance(usecols, list):
        usecols = [usecols]

    # close connection
    cursor.close()

    return pd.DataFrame(data, columns=usecols)
