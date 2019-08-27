#!/usr/bin/python3

'''
Prerequisites: 
	- the header of every sequence in the MSA file must start with a label that can identify every FunFam group and metagenome sample group,
		and an underscore must immediately follow the label i.e. >[label]_[seq_id]
	- the description line of the fasta sequence should contain the FunFam id at the end in the format [C].[A].[T].[H]/FF/[FF id]
# output: screen output of individual list of mutant residues in each FunFam
# also output a table of MG seqs with novel mutations
'''

import sys
import re
import msa_position_converter as mpc
import pandas as pd
from Bio import AlignIO
from collections import Counter

# Command-line arguments
msa_file = sys.argv[1]  # fasta MSA with all funfam seqs and metagenome seqs together
ref_seq_id = sys.argv[2]  # must be the id of a sequence in the msa file
funfam_group_ids = sys.argv[3]  # comma-separated list (w/o spaces) of one or more strings that are at the start of every known funfam sequence
novel_group_ids = sys.argv[4]  # comma-separated list (w/o spaces) of one or more strings that are at the start of every novel/uncharacterised sequence
ref_sdp_positions = sys.argv[5]  # comma-separated list (w/o spaces) of one or more e.g. 96,142,248,257

# Reads in the multiple sequence alignment as a Bio.AlignIO object
msa_object = AlignIO.read(msa_file, 'fasta')

# Allows a shorter version of the ref_seq id to find the full version of the id (which includes start position of domain seq)
for record in msa_object:
	if ref_seq_id in record.id:
		ref_seq_fullid = record.id
# Creates dictionary of corresponding positions with ref_seq positions as keys and MSA positions as values
ref_to_msa_dict = mpc.ref_to_msa(msa_object, ref_seq_id, start=int(ref_seq_fullid.split('/')[1].split('-')[0]))

# Creates list of sdp positions in the ref seq based on the user input (they will be translated to msa positions when needed)
ref_sdp_positions = [int(x) for x in ref_sdp_positions.split(',')]


funfam_group_ids = funfam_group_ids.split(',')
novel_group_ids = novel_group_ids.split(',')
p_group = re.compile(r'^(.+?)_')

# Initialise pandas dataframe with row names as sequence ids and column names as sdp positions (ref seq).
# First row will be called 'FF consensus' and have a list of all AAs at that position in FunFam sequences.
# It will be filled with sequences that have at least one novel mutation.
mutant_df = pd.DataFrame(columns=sorted(ref_sdp_positions))
mutant_df_summary = pd.DataFrame(columns=sorted(ref_sdp_positions)) # will only be filled with two rows


for position in ref_sdp_positions:
	ff_groups_consensus = []
	for record in msa_object:
		group = p_group.search(record.id).group(1)
		if group in funfam_group_ids:
			ff_groups_consensus.append(record.seq[ref_to_msa_dict[position] - 1])

	ff_groups_consensus = set(ff_groups_consensus)
	mutant_df.loc['FF consensus', position] = ''.join(ff_groups_consensus)
	mutant_df_summary.loc['FF consensus', position] = ''.join(ff_groups_consensus)
	mutant_df_summary.loc['Mutations', position] = ''

	# 2nd iteration over the sequences, this time to find novel mutants
	for record in msa_object:
		group = p_group.search(record.id).group(1)
		if group in novel_group_ids:
			residue = record.seq[ref_to_msa_dict[position] - 1]
			if (residue not in ff_groups_consensus) and (residue != '-'):
				mutant_df.loc[record.id, position] = residue
				mutant_df_summary.loc['Mutations', position] += residue

for position in mutant_df:
	mutant_df.loc['FF consensus', position] = ''.join(sorted(mutant_df.loc['FF consensus', position]))

# Mutations sorted in order of frequency
for position in mutant_df_summary:
	counts = Counter(list(mutant_df_summary.loc['Mutations', position]))
	mutant_df_summary.loc['Mutations', position] = ''.join(sorted(list(mutant_df_summary.loc['Mutations', position]), key=lambda x: -counts[x]))

print(mutant_df.to_csv(sep='\t'))  # convert to table in a word document - the output on the screen does not align the columns correctly
print(mutant_df_summary.to_csv(sep='\t'))


#############################################################################
## TO DO

# the novel mutations are not separated by closest funfams
# separate tables can be made by using each funfam separately and only the weak hits that most closely match that funfam
# compare e-value score for weak hits with and without novel mutations
# calculate co-variation

# watch out for short fragments that aren't aligned well (delete them in jalview beforehand)
