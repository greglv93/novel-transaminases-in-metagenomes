#!/usr/bin/python3

'''
Takes in crh files (genomescan output) for one or more samples.
Generates a series of reports (stats tables) on FunFam hits within a specified CATH Superfamily.
'''

import re 
import operator
import argparse


# Creates parser object and defines optional arguments
parser = argparse.ArgumentParser(description='Takes in crh files and generates a series of reports', add_help=True)
parser.add_argument('-s', '--samples', dest='samples', required=True, nargs='+', help='Enter one or more sample names')
parser.add_argument('-f', '--supfams', dest='supfams', required=True, nargs='+', help='Enter one or more CATH Superfamilies e.g. 3.40.640.10')
# Reads in command-line argument and adds the supfam attribute to the parser object
args = parser.parse_args()


# Regular expressions to extract information from each line of the crh file
p_gene_id = re.compile(r'^(\w+?[+-])')
p_cath = re.compile(r'^\S+?\s(.+?)\/')
p_ff = re.compile(r'^\S+?\s\S+?\/FF\/(\d+)')
p_gene_bounds = re.compile(r'^k141_\d+?_(\d+?)_(\d+)_')
p_dom_bounds = re.compile(r'^(?:\S+?\s){4}(\d+?)-(\d+?)\s')
p_eval = re.compile(r'^(?:\S+?\s){5}(.+?)\s')  # Conditional E-value
p_score = re.compile(r'^(?:\S+?\s){2}(.+?)\s')


###############################################
print('reading in inclusion thresholds...')
###############################################

inc_thresh = {}
for supfam in args.supfams:

	with open('./inc_thresh/' + supfam, 'r') as f:
		data = f.read().splitlines()
	for i in range(0, len(data)):
		data[i] = data[i].split('\t')
		data[i][0] = supfam + '/' + data[i][0]
	# index 0 = funfam id; index 2 = e-value threshold
	data_dict = {k[0]:k[2] for k in data}
	inc_thresh.update(data_dict)

'''
inc_thresh_score = {}
for supfam in args.supfams:

	with open(inc_thresh_path + supfam, 'r') as f:
		data = f.read().splitlines()
	for i in range(0, len(data)):
		data[i] = data[i].split('\t')
		data[i][0] = supfam + '/' + data[i][0]
	# index 0 = funfam id; index 1 = score threshold
	data_dict = {k[0]:k[1] for k in data}
	inc_thresh_score.update(data_dict)
'''

###############################################

ff_count_strong_fullgenes_allsamples = {}
ff_count_strong_genefrags_allsamples = {}
ff_count_weak_fullgenes_allsamples = {}
ff_count_weak_genefrags_allsamples = {}

###############################################

for sample in args.samples:

	###############################################
	print('\nreading in %s crh file...' % sample)
	###############################################

	with open('../genomescan/results/sample/run2/sample_genes_run1_hits.crh'.replace('sample', sample), 'r') as f:
		crh = f.read().splitlines()
	crh = list(filter(None, crh[2:]))  # removes the header and any blank lines in the crh file
	total_hits = len(crh)  # will be used to measure progress

	###############################################
	print('reading in full-length gene IDs...')
	###############################################

	with open('./' + sample + '_full_gene_ids', 'r') as f:
		all_full_genes = f.read().splitlines()

	for i in range(0, len(all_full_genes)):
		all_full_genes[i] = all_full_genes[i].lstrip('>')


	###############################################
	print('\ncounting all %s hits: ' % sample)
	###############################################

	ff_count_strong_fullgenes = {}
	ff_count_strong_genefrags = {}
	ff_count_weak_fullgenes = {}
	ff_count_weak_genefrags = {}

	strong_hits_fullgenes = []
	strong_hits_genefrags = []
	weak_hits_fullgenes = []
	weak_hits_genefrags = []
	domainfrags = []

	for hit in crh:

		gene = p_gene_id.search(hit).group(1)
		cath = p_cath.search(hit).group(1)
		ff = p_ff.search(hit).group(1)
		domain_id = cath + '/' + ff

		if cath in args.supfams:

			if gene in all_full_genes:
				full_gene = True
			else:
				full_gene = False
				gene_bounds = p_gene_bounds.search(hit)
				gene_len = (int(gene_bounds.group(2)) - int(gene_bounds.group(1))) + 1
				dom_bounds = p_dom_bounds.search(hit)
				dom_len = (int(dom_bounds.group(2)) - int(dom_bounds.group(1))) + 1
				if (int(dom_bounds.group(1)) <= 5) or (gene_len - int(dom_bounds.group(2)) < 5) or (dom_len < 200):
					domainfrags.append(hit)
					# not foolproof - a domain could be incomplete but the hmm decides the best matching window doesn't start right at the boundary
					continue

			# Strong hits are counted when the E-value meets the inclusion threshold of the FunFam profile
			if float(p_eval.search(hit).group(1)) <= float(inc_thresh[domain_id]):
				if full_gene:
					strong_hits_fullgenes.append(hit)
					if domain_id not in ff_count_strong_fullgenes:
						ff_count_strong_fullgenes[domain_id] = 1
					else:
						ff_count_strong_fullgenes[domain_id] += 1
					if domain_id not in ff_count_strong_fullgenes_allsamples:
						ff_count_strong_fullgenes_allsamples[domain_id] = 1
					else:
						ff_count_strong_fullgenes_allsamples[domain_id] += 1
				elif not full_gene:
					strong_hits_genefrags.append(hit)
					if domain_id not in ff_count_strong_genefrags:
						ff_count_strong_genefrags[domain_id] = 1
					else:
						ff_count_strong_genefrags[domain_id] += 1
					if domain_id not in ff_count_strong_genefrags_allsamples:
						ff_count_strong_genefrags_allsamples[domain_id] = 1
					else:
						ff_count_strong_genefrags_allsamples[domain_id] += 1


			# Weak hits are counted when the E-value does not meet the inclusion threshold of the FunFam profile
			elif 0.001 >= float(p_eval.search(hit).group(1)) > float(inc_thresh[domain_id]):
				# This represents domains that do not strongly match any known FunFam,
				# but the closest match is a transaminase funfam;
				# These could represent new funfams or slightly different substrate specificity within an existing funfam
				if full_gene:
					weak_hits_fullgenes.append(hit)
					if domain_id not in ff_count_weak_fullgenes:
						ff_count_weak_fullgenes[domain_id] = 1
					else:
						ff_count_weak_fullgenes[domain_id] += 1
					if domain_id not in ff_count_weak_fullgenes_allsamples:
						ff_count_weak_fullgenes_allsamples[domain_id] = 1
					else:
						ff_count_weak_fullgenes_allsamples[domain_id] += 1
				elif not full_gene:
					weak_hits_genefrags.append(hit)
					if domain_id not in ff_count_weak_genefrags:
						ff_count_weak_genefrags[domain_id] = 1
					else:
						ff_count_weak_genefrags[domain_id] += 1
					if domain_id not in ff_count_weak_genefrags_allsamples:
						ff_count_weak_genefrags_allsamples[domain_id] = 1
					else:
						ff_count_weak_genefrags_allsamples[domain_id] += 1

		'''
		if cath in args.supfams:

			if float(p_score.search(hit).group(1)) >= float(inc_thresh_score[domain_id]):

				if domain_id not in ff_count_strong:
					ff_count_strong[domain_id] = 1
				else:
					ff_count_strong[domain_id] += 1


			elif float(p_score.search(hit).group(1)) < float(inc_thresh_score[domain_id]):
				# this represents domains that do not strongly match any known funfam,
				# but the closest match is a transaminase funfam;
				# these could represent new funfams or slightly different substrate specificity within an existing funfam

				if domain_id not in ff_count_weak:
					ff_count_weak[domain_id] = 1
				else:
					ff_count_weak[domain_id] += 1
		'''

	ff_count_strong = { k: ff_count_strong_fullgenes.get(k, 0) + ff_count_strong_genefrags.get(k, 0) \
						for k in set(ff_count_strong_fullgenes) | set(ff_count_strong_genefrags) }
	ff_count_weak = { k: ff_count_weak_fullgenes.get(k, 0) + ff_count_weak_genefrags.get(k, 0) \
						for k in set(ff_count_weak_fullgenes) | set(ff_count_weak_genefrags) }


	###############################################
	# stats tables
	###############################################

	dicts_list = [ [ff_count_strong_fullgenes, 'strong', 'fullgenes'], [ff_count_strong_genefrags, 'strong', 'genefrags'], 
					[ff_count_strong, 'strong', 'all'], [ff_count_weak_fullgenes, 'weak', 'fullgenes'], 
					[ff_count_weak_genefrags, 'weak', 'genefrags'], [ff_count_weak, 'weak', 'all'] ]

	for ff_count_dict in dicts_list:

		sorted_ff_count_dict = sorted(ff_count_dict[0].items(), key=operator.itemgetter(1), reverse=True)

		with open('./output/supfam/sample_stats_strength'.replace('supfam', \
				'&'.join(args.supfams)).replace('sample', sample).replace('strength', ff_count_dict[1]) \
				+ ff_count_dict[2], 'w') as out:
			header = 'SF\tFF\thits\n'
			out.write(header)
			for ff in sorted_ff_count_dict:
				line = ff[0].split('/')[0] + '\t' + ff[0].split('/')[1] + '\t' + str(ff[1]) + '\n'
				out.write(line)

	###############################################
	# filtered crh files
	###############################################

	crh_lists = [[strong_hits_fullgenes, 'strong', 'fullgenes'], [strong_hits_genefrags, 'strong', 'genefrags'], \
					[weak_hits_fullgenes, 'weak', 'fullgenes'], [weak_hits_genefrags, 'weak', 'genefrags']]

	for crh_list in crh_lists:

		with open('./output/supfam/sample_hits_strength'.replace('supfam', \
					'&'.join(args.supfams)).replace('sample', sample).replace('strength', crh_list[1]) \
					+ crh_list[2] + '.crh', 'w') as out:
			for hit in crh_list[0]:
				out.write(hit + '\n')

	###############################################
	# domain fragments list
	###############################################

	with open('./output/supfam/sample_domainfrags.crh'.replace('supfam', \
					'&'.join(args.supfams)).replace('sample', sample), 'w') as out:
		for hit in domainfrags:
			out.write(hit + '\n')
	###############################################


###################################################
# combined stats tables for all samples
###################################################

ff_count_strong_allsamples = { k: ff_count_strong_fullgenes_allsamples.get(k, 0) + ff_count_strong_genefrags_allsamples.get(k, 0) \
					for k in set(ff_count_strong_fullgenes_allsamples) | set(ff_count_strong_genefrags_allsamples) }
ff_count_weak_allsamples = { k: ff_count_weak_fullgenes_allsamples.get(k, 0) + ff_count_weak_genefrags_allsamples.get(k, 0) \
					for k in set(ff_count_weak_fullgenes_allsamples) | set(ff_count_weak_genefrags_allsamples) }

dicts_list = [ [ff_count_strong_fullgenes_allsamples, 'strong', 'fullgenes'], [ff_count_strong_genefrags_allsamples, 'strong', 'genefrags'], \
				[ff_count_strong_allsamples, 'strong', 'all'], [ff_count_weak_fullgenes_allsamples, 'weak', 'fullgenes'], \
				[ff_count_weak_genefrags_allsamples, 'weak', 'genefrags'], [ff_count_weak_allsamples, 'weak', 'all'] ]

for ff_count_dict in dicts_list:

	sorted_ff_count_dict = sorted(ff_count_dict[0].items(), key=operator.itemgetter(1), reverse=True)

	with open('./output/supfam/allsamples_stats_strength'.replace('supfam', \
			'&'.join(args.supfams)).replace('strength', ff_count_dict[1]) \
			+ ff_count_dict[2], 'w') as out:
		header = 'SF\tFF\thits\n'
		out.write(header)
		for ff in sorted_ff_count_dict:
			line = ff[0].split('/')[0] + '\t' + ff[0].split('/')[1] + '\t' + str(ff[1]) + '\n'
			out.write(line)
