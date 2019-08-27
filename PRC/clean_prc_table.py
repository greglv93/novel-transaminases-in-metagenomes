#!/usr/bin/python3

'''
Prerequisites:
    - output from iterative_PRC.sh (concatenation of all of the 'one vs all' PRC runs).
        (This will contain redundant headers and 'END' statements that will be cleaned here)
    - SC_FF mapping file (in ../mappings/), including columns for FF size and modal EC annotation
    - run ../Scripts/create_funfam_table.py with no optional arguments to provide table that can be merged with SC_FF table
Does some cleaning and wrangling to produce a simple network table and a more complicated one with several node attributes.
The latter can be imported into Cytoscape for network visualisation.
'''

import argparse
import re
import pandas as pd

# Creates parser object and defines optional arguments
parser = argparse.ArgumentParser()
parser.add_argument('-s', '--supfam', dest='supfam', required=True, help='give superfamily of interest e.g. 3.40.640.10')
# Reads in command-line argument and adds the supfam attribute to the parser object
args = parser.parse_args()

# Pattern for extracting the 3 columns of interest
p = re.compile(r'(\S+)\s+(?:\S+\s+){4}(\S+)(?:\s+\S+){5}\s+(\S+)')
# group(1) = FF1
# group(2) = FF2
# group(3) = reverse score

# Reads in messy concatenation of PRC output and extracts a list of the rows of data
with open('./global_global_PRC', 'r') as f:
    lines = f.read().splitlines()
header = lines[13] + '\n'
data = []
for l in lines:
    if args.supfam in l:  # this ensures only rows with data are included
        data.append(l)

# Writes full cleaned table
with open('./gg_prc_clean', 'w') as out:
    out.write(header)
    out.write('\n'.join(data))

# Removes duplicates (FF1 vs FF2 and FF2 vs FF1) and writes basic table with only the two FunFams and their score
with open('./gg_prc_clean_noduplicates', 'w') as out:
    out.write('hmm1\thmm2\talignment_score\n')
    pairs = {}
    for x in data:
        ff1 = p.search(x).group(1).split('/FF/')[1]
        ff2 = p.search(x).group(2).split('/FF/')[1]
        score = p.search(x).group(3)

        if ff1 not in pairs:
            if ff2 not in pairs:
                pairs[ff1] = [ff2]
                out.write(ff1 + '\t' + ff2 + '\t' + score + '\n')
            else:
                if ff1 not in pairs[ff2]:
                    pairs[ff2].append(ff1)
                    out.write(ff1 + '\t' + ff2 + '\t' + score + '\n')
                else:
                    pass
        else:
            if ff2 not in pairs[ff1]:
                pairs[ff1].append(ff2)
                out.write(ff1 + '\t' + ff2 + '\t' + score + '\n')
            else:
                pass

# Merges SC_FF mapping table with FunFam table
sc_table = pd.read_csv('../mappings/' + args.supfam + '.SC_FF.mapping_withEC3and4mode.txt', sep='\t')
sc_table = sc_table.set_index('FF')
sc_table = sc_table.rename(columns={'SSG5_v4.2':'SC5'})
sc_table = sc_table.drop(['SF', 'SSG9_v4.1', 'REP', 'FF_SIZE', 'modal_EC4', 'modal_EC3'], axis=1)
full_table = pd.read_csv('./' + args.supfam + '-v4_1_0-funfams-_uniprot', sep='|', skiprows=2,
    names=['FF','FF_size','DOPS','REP','targetEC3perc','EC4_mode', 'EC4_mode_perc', 'FF_name'])
full_table = full_table.set_index('FF')
full_table = full_table.drop(['REP','targetEC3perc','EC4_mode_perc','FF_name'], axis=1)
node_table = pd.merge(sc_table, full_table, how='outer', left_index=True, right_index=True)
EC4_mode = list(node_table.loc[:,'EC4_mode'])
EC3_mode = ['.'.join(x.split('.')[:-1]) for x in EC4_mode]
node_table['EC3_mode'] = pd.Series(EC3_mode, index=node_table.index)

# Reads in scores from basic table (created 2 blocks above)
with open('../PRC/gg_prc_clean_noduplicates', 'r') as f:
    f_lines = f.read().splitlines()[1:]
prc_scores = [x.split('\t') for x in f_lines]

# Writes full network table with several attributes for each FunFam
with open('../PRC/3.40.640.10_FFs_prc_network', 'w') as out:
    out.write('FF1\tSC5_1\tFF_size_1\tmodal_EC3_1\tFF2\tSC5_2\tFF_size_2\tmodal_EC3_2\tprc_score\n')
    for edge in prc_scores:
        out.write(edge[0] + '\t' + 
                    str(node_table.loc[float(edge[0]),'SC5']) + '\t' + 
                    str(node_table.loc[float(edge[0]), 'FF_size'])  + '\t' + 
                    str(node_table.loc[float(edge[0]), 'EC3_mode']) + '\t' + 
                    edge[1] + '\t' + 
                    str(node_table.loc[float(edge[1]),'SC5']) + '\t' + 
                    str(node_table.loc[float(edge[1]), 'FF_size'])  + '\t' + 
                    str(node_table.loc[float(edge[1]), 'EC3_mode']) + '\t' +
                    edge[2] + '\n'
                    )
