'''
Module that can be used to translate between sequence positions in a multiple sequence alignment and a reference sequence within the alignment.
Used in sort_groupsim_scores.py and novel_mutant_finder.py
Example situation where this module is required: 
	- a tool or function returns data associated with positions in MSA1;
	- further sequences are then added to MSA1 to create MSA2
	- analysis needs to be done on MSA2 based on MSA1 positions
	- a reference sequence, present in both MSAs, can be used to bridge the MSAs by providing corresponding positions
'''

def msa_to_ref(msa_object, ref_seq_id, start=1):
	# Returns dictionary of corresponding positions with msa positions as keys and ref_seq positions as values (gap values are strings).
	# msa_object must be a biopython sequene object (use Bio.AlignIO.read(file, 'fasta') in the scripts where this function is used)
		# class 'Bio.Align.MultipleSeqAlignment'
	# Note: ref_seq just means reference sequence, it does not refer to the NCBI database
	# start is the start position of the sequence 
		# change to the start of a domain or first position of a truncated sequence if you want the position on the full sequence

	for record in msa_object:
		if ref_seq_id in record.id:  # the ref_seq_id argument can be an incomplete id, as long as it is unique
			ref_seq = str(record.seq)
			break

	position_dict = {}
	current_msa_position = 1
	current_ref_position = start - 1
	while current_msa_position <= len(ref_seq):
		if ref_seq[current_msa_position - 1] == '-':
			position_dict[current_msa_position] = 'gap:' + str(current_ref_position) + '-' + str(current_ref_position + 1)
		else:
			current_ref_position += 1
			position_dict[current_msa_position] = current_ref_position
		current_msa_position += 1

	return position_dict


def ref_to_msa(msa_object, ref_seq_id, start=1):
	# Returns dictionary of corresponding positions with ref_seq positions as keys and msa positions as values (all int)

	for record in msa_object:
		if ref_seq_id in record.id:
			ref_seq = str(record.seq)
			break

	position_dict = {}
	current_msa_position = 1
	current_ref_position = start - 1
	while current_msa_position <= len(ref_seq):
		if ref_seq[current_msa_position - 1] != '-':
			current_ref_position += 1
			position_dict[current_ref_position] = current_msa_position
		else:
			pass
		current_msa_position += 1

	return position_dict

# Both functions were tested, using the output of one as the input of the other (translating and translating back); and manual check of an alignment viewer
	# TO DO: prove this with a unit test

# A copy of this module was made in /usr/local/lib/python3.5/dist-packages so it can be imported from anywhere


## TO DO:
# change comments to docstrings
	# mention that the functions return a list of translated positions in the order of the input positions
# make the functions return an error if the record id is not found in the msa object
