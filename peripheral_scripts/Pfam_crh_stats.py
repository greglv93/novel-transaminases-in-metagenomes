#!/usr/bin/python3

'''
For a given set of .crh files, this script generates 4 reports: CATH superfamily composition, funfam composition, 
CATH superfamily architectures, and funfam architectures
'''

import re
import operator

PFAM_IDS = ['PF00155', 'PF00202', 'PF01063', 'PF00266', 'PF01041', 'PF12897']

# regex patterns for important parts of the '.crh' data format
p_uniprotid = re.compile(r'^(\w+?)\s')  # captures Uniprot ID
p_cath = re.compile(r'^\S+?\s(.+?)\/')  # captures the cath id
p_ff = re.compile(r'^\S+?\s\S+?\/(FF\/\d+)')  # captures the funfam id


ff_count_all_pfams = {}

for pfamid in PFAM_IDS:

	with open('./' + pfamid + '_full_length_sequences.crh', 'r') as file:
		crh = file.read().splitlines()
	
	# removes blank strings
	crh = list(filter(None, crh[2:]))

	# simple count of cath superfamily composition of the .crh file
	cath_count = {}
	for hit in crh:
		if p_cath.search(hit) is not None:
			cath_id = p_cath.search(hit).group(1)
			if cath_id not in cath_count:
				cath_count[cath_id] = 1
			else:
				cath_count[cath_id] += 1
		else:
			pass
		
	sorted_cath_count = sorted(cath_count.items(), key=operator.itemgetter(1), reverse=True)

	with open('./' + pfamid + '_crh_cathcomposition', 'w') as out:
		for x in sorted_cath_count:
				out.write(str(x[1]) + '\t' + x[0] + '\n')


	# simple count of CATALYTIC funfam composition of the .crh file
		# also adds up the no. of genes with the desired domain; subtracting from total will give the structurally impure section of the pfam
	ff_count = {}
	uniprot_ids = []  # only those with a catalytic TAm domain
	all_uniprot_ids = []
	for x in crh:
		if p_cath.search(x) is not None:
			all_uniprot_ids.append(p_uniprotid.search(x).group(1))
			if (p_cath.search(x).group(1) == '3.20.10.10') or (p_cath.search(x).group(1) == '3.40.640.10'):
				full_id = p_cath.search(x).group(1) + '.' + p_ff.search(x).group(1)
				if full_id not in ff_count:
					ff_count[full_id] = 1
					ff_count_all_pfams[full_id] = 1
				else:
					ff_count[full_id] += 1
					ff_count_all_pfams[full_id] += 1
				uniprot_id = p_uniprotid.search(x).group(1)
				if uniprot_id not in uniprot_ids:
					uniprot_ids.append(uniprot_id)
		else:
			pass

	sorted_ff_count = sorted(ff_count.items(), key=operator.itemgetter(1), reverse=True)

	with open('./' + pfamid + '_crh_catalytic_funfam_composition', 'w') as out:
		for x in sorted_ff_count:
			out.write(str(x[1]) + '\t' + x[0] + '\n')

	all_uniprot_ids = set(all_uniprot_ids)
	print('%s seqs with a catalytic TAm domain: %d (%.2f percent)' % (pfamid, len(uniprot_ids), 
		(float(len(uniprot_ids))/float(len(all_uniprot_ids)))*100)
		)


	# creates lists of architectures in format [[uniprot_id, domain_1, domain_2], [etc]]
	last_uniprot_id = p_uniprotid.search(crh[0]).group(1)
	cath_architectures = [[last_uniprot_id]]
	funfam_architectures = [[last_uniprot_id]]
	for x in crh:
		if p_uniprotid.search(x) is not None:
			current_uniprot_id = p_uniprotid.search(x).group(1)
			cath_id = p_cath.search(x).group(1)
			ff_id = p_ff.search(x).group(1).replace('/', '')
		else:
			continue

		if current_uniprot_id == last_uniprot_id:
			cath_architectures[-1].append(cath_id)
			funfam_architectures[-1].append(cath_id + '.' + ff_id)
		else:
			cath_architectures.append([current_uniprot_id])
			cath_architectures[-1].append(cath_id)
			funfam_architectures.append([current_uniprot_id])
			funfam_architectures[-1].append(cath_id + '.' + ff_id)

		last_uniprot_id = current_uniprot_id

	# generates dicts of architecture counts and writes reports to output files	
	cath_architecture_count = {}
	for x in cath_architectures:
		if '---'.join(x[1:]) not in cath_architecture_count:
			cath_architecture_count['---'.join(x[1:])] = 1
		else:
			cath_architecture_count['---'.join(x[1:])] += 1
	sorted_cath_architecture_count = sorted(cath_architecture_count.items(), key=operator.itemgetter(1), reverse=True)
	sorted_cath_architecture_count_table = ''
	for x in sorted_cath_architecture_count:
		sorted_cath_architecture_count_table += (str(x[1]) + '\t' + str(x[0]) + '\n')

	with open('./' + pfamid + '_crh_cath_architectures', 'w') as out:
		out.write('\n' + '-'*40 + '\narchitecture composition:\n' + '-'*40 + '\n' + sorted_cath_architecture_count_table)
		out.write('\n' + '-'*40 + '\nUniprot ID\tarchitecture\n' + '-'*40 + '\n')
		for x in cath_architectures:
			out.write(x[0] + '\t' + '---'.join(x[1:]) + '\n')
						

	funfam_architecture_count = {}
	for x in funfam_architectures:
		if '---'.join(x[1:]) not in funfam_architecture_count:
			funfam_architecture_count['---'.join(x[1:])] = 1
		else:
			funfam_architecture_count['---'.join(x[1:])] += 1
	sorted_funfam_architecture_count = sorted(funfam_architecture_count.items(), key=operator.itemgetter(1), reverse=True)
	sorted_funfam_architecture_count_table = ''
	for x in sorted_funfam_architecture_count:
		sorted_funfam_architecture_count_table += (str(x[1]) + '\t' + str(x[0]) + '\n')

	with open('./' + pfamid + '_crh_funfam_architectures', 'w') as out:
		out.write('\n' + '-'*40 + '\narchitecture composition:\n' + '-'*40 + '\n' + sorted_funfam_architecture_count_table)
		out.write('\n' + '-'*40 + '\nUniprot ID\tarchitecture\n' + '-'*40 + '\n')
		for x in funfam_architectures:
			out.write(x[0] + '\t' + '---'.join(x[1:]) + '\n')


# count of CATALYTIC funfam composition of all transaminase Pfams
sorted_ff_count_all_pfams = sorted(ff_count_all_pfams.items(), key=operator.itemgetter(1), reverse=True)

with open('./allPfams_crh_catalytic_funfam_composition', 'w') as out:
	for x in sorted_ff_count_all_pfams:
		out.write(str(x[1]) + '\t' + str(round((x[1]/sum(ff_count_all_pfams.values()))*100, 2)) + '%' + '\t' + x[0] + '\n')
