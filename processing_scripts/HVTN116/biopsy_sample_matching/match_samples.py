## ---------------------------------------------------------------------------##
# Author: Sara Thiebaud
# Date: 9 June 2026
# Purpose: HVTN 116 sample metadata matching for the McElrath lab, request by Yunda and Lily
## ---------------------------------------------------------------------------##
import pandas
import numpy as np
import os
import datetime as datetime
import pdb

import sdmc_tools.process as sdmc_tools
import sdmc_tools.access_ldms as access_ldms

NETWORK = 'hvtn'
PROTOCOL = 116

WORKING_DIR = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN116/assays/biopsy_sample_matching/'
LAB_METADATA_FILE = 'HVTN116 MUC BIOPSY META DATA FROM BINDERS BATCH 12-64 NW 2.xlsx'
QDATA_FOLDER = '/trials/vaccine/p116/s001/qdata/'
QDATA_FILE = 'VTN116_Biopsy Mass_20230208.txt'

ldms = access_ldms.pull_one_protocol(NETWORK, PROTOCOL)
ldms['guspec_core'] = ldms.guspec.str.rpartition("-")[0]
ldms = ldms.drop(columns='guspec').drop_duplicates()

biopsy_mass = pandas.read_csv(os.path.join(QDATA_FOLDER, QDATA_FILE), dtype=str, sep='\t')
biopsy_mass['guspec_core_lysate'] = biopsy_mass.apply(lambda row: row['guspec'].rsplit('-', 1)[0], axis=1)
biopsy_mass['SAMPLE TYPE'] = biopsy_mass.apply(lambda row: 'VAG' if 'Vaginal' in row['spectype'] else 'CER' if 'Cervical' in row['spectype'] else '', axis=1)

lab_metadata = pandas.read_excel(os.path.join(WORKING_DIR, LAB_METADATA_FILE))
lab_metadata.dropna(subset=['Lysate Spec IDs'], inplace=True)
lab_metadata['guspec_core_lysate'] = lab_metadata.apply(lambda row: row['Lysate Spec IDs'].rsplit('-', 2)[0] if '(' in row['Lysate Spec IDs'] else row['Lysate Spec IDs'].rsplit('-', 1)[0], axis=1)
pdb.set_trace()
lab_ds = pandas.merge(lab_metadata, ldms, how='left', left_on='guspec_core_lysate', right_on='guspec_core', indicator=True)
pdb.set_trace()

#m = pandas.merge(biopsy_mass, lab_ds, how='outer', left_on='guspec_core_lysate', right_on=[', indicator=True)
m = pandas.merge(biopsy_mass, lab_ds, how='outer', left_on=['ptid', 'visitno', 'SAMPLE TYPE'], right_on=['txtpid', 'vidval', 'SAMPLE TYPE'], indicator=True)
pdb.set_trace()
m.to_csv('HVTN116_biopsy_mass_with_lysate_volume_added_by_sdmc_2026-06-09.txt', index=False, sep='\t')
pdb.set_trace()