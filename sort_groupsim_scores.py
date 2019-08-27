#!/usr/bin/python3

'''
Sorts GroupSim output and writes top scores (user-defined percentage) and their positions - in both the MSA and a given reference sequence - to a new file
Also writes a Chimera select command (to the end of the output file) for the top scoring residues
USAGE: ./sort_groupsim_scores.py [multiple sequence alignment] [reference sequence id] [groupsim output] [percentage] [output file]
Example: ./sort_groupsim_scores.py ./MSA1.fasta 3du4A02 ./groupsim_scores_MSA1 10 ./groupsim_scores_MSA1_top10pc
Important: the multiple sequence alignment file must be the same input as was given to group_sim_sdp.py,
			otherwise the translation of sequence positions may be wrong
'''

import sys
import re
import operator
import msa_position_converter as mpc
from Bio import AlignIO

# Command-line arguments
msa_file = sys.argv[1]  # must be in fasta format
ref_seq_id = sys.argv[2]
groupsim_output = sys.argv[3]
percentage = int(sys.argv[4])  # must be integer value
output_file = sys.argv[5]

# Reads multpiple sequence alignment file
msa_object = AlignIO.read(msa_file, 'fasta')
# Allows a shorter version of the ref_seq id to find the full version of the id (which includes start position of domain seq)
for record in msa_object:
	if ref_seq_id in record.id:
		ref_seq_fullid = record.id
# Creates dictionary of corresponding positions with msa positions as keys and ref_seq positions as values
msa_to_ref_dict = mpc.msa_to_ref(msa_object, ref_seq_id, start=int(ref_seq_fullid.split('/')[1].split('-')[0]))

# Reads groupsim output file, removing the header lines and any blank lines
with open(groupsim_output, 'r') as f:
	msa_cols = f.read().splitlines()[4:]
msa_cols = list(filter(None, msa_cols))

# Patterns for extracting position and score from each line of the groupsim output
p_position = re.compile(r'^(\d+)\t')
p_score = re.compile(r'^\d+\t(\S+)\t')

# Puts position and score of each column in the MSA into a dict
scores_dict = {}
for x in msa_cols:
	msa_position = int(p_position.search(x).group(1)) + 1  # groupsim starts from position 0
	score = p_score.search(x).group(1)
	if score == 'None':
		score = 0
	else:
		score = round(float(score), 5)
	scores_dict[msa_position] = score

# Sorts the dict by score, in descending order
sorted_scores_dict = sorted(scores_dict.items(), key=operator.itemgetter(1), reverse=True)

# Calculates the average sequence length of the MSA
total_seqs = len(''.join(msa_cols[0].split('\t')[-1].split(' | ')))
total_nongap_chars = 0
for x in msa_cols:
	chars = ''.join(x.split('\t')[-1].split(' | '))
	nongap_chars = chars.replace('-', '')
	total_nongap_chars += len(nongap_chars)
avg_seq_len = total_nongap_chars / total_seqs

# No. of top scores to display based on user input
scores_to_display = round((avg_seq_len)*(percentage/100))

# Writes output file
with open(output_file, 'w') as out:
	out.write('# input file: %s \n' % groupsim_output)
	out.write('# average sequence length in MSA: %d\n' % avg_seq_len)
	out.write('# top %d percent = top %s groupsim scores:\n' %(percentage, scores_to_display))
	out.write('Score\tMSA_col\t%s_position\t(note: MSA_col starts from 1, while groupsim col_num starts from 0)\n' % ref_seq_id)
	for x in sorted_scores_dict[:scores_to_display]:
		out.write(str(x[1]) + '\t' + str(x[0]) + '\t' + str(msa_to_ref_dict[x[0]]) + '\n')
	# Writes chimera select command
	out.write('\nChimera select command:\nsel :')
	for x in sorted_scores_dict[:scores_to_display-1]:
		out.write(str(msa_to_ref_dict[x[0]]) + ',')
	out.write(str(msa_to_ref_dict[sorted_scores_dict[scores_to_display-1][0]]) + '\n')
