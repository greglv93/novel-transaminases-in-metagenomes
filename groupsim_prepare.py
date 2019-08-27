#!/usr/bin/python3

'''
Reads a fasta file and appends a group name to the end of each sequence, to meet the required format for GroupSim.
The group name is also added to the start of each fasta header, to allow easier manual inspection of several groups of sequences in jalview.

USAGE:
./groupsim_prepare.py [input fasta] [group name] [output file]

'''

import sys
from Bio import SeqIO

# Command-line arguments
input_fasta = sys.argv[1]
group_name = sys.argv[2]
output_fasta = sys.argv[3]

# appends group name to start of headers for visibility and end of headers to meet the required format for GroupSim
with open(output_fasta, 'w') as out:
	for record in SeqIO.parse(input_fasta, 'fasta'):
		record.description = group_name + '_' + record.description + ' |' + group_name
		out.write('>' + str(record.description) + '\n')
		out.write(str(record.seq) + '\n')
