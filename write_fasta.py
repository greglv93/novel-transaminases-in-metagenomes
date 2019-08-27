#!/usr/bin/python3

'''
Takes in a large fasta file (superset of sequences of interest) and crh file with hits of interest (produced by genomescan_stats.py).
Takes in command-line arguments that specify which crh file to use and which FunFam to look for.
Uses crh file to extract the domain sequences (not full gene sequences) of interest.
Appends the FunFam ID, domain boundaries, and E-value to the fasta headers.
'''

import argparse
import sys
import re
import itertools

# Creates parser object and defines optional arguments
parser = argparse.ArgumentParser(description='Writes fasta amino acid files for FunFam hits in metagenomic sequences, \
									with options to include only strong or weak hits and to use only full genes or fragments', add_help=True)
parser.add_argument('-s', '--samples', dest='samples', required=True, nargs='+', help='list of one or more samples e.g. "pilluana maras3"')
parser.add_argument('-j', '--joined', action='store_true', help='if -j option is given, the sample sequences will be joined into one output fasta, \
					else, separate output fastas will be created for each sample (if more than one given)')
parser.add_argument('-f', '--funfam', dest='funfam', required=True, help='full FunFam code (e.g. 3.40.640.10/FF/63207)')
parser.add_argument('-t', '--type', dest='type', default='all', help='"strong", "weak", or "all"')
parser.add_argument('-l', '--length', dest='length', default='fullgenes', help='"fullgenes" or "genefrags"')

# Reads in command-line arguments and adds the attributes to the parser object
args = parser.parse_args()

# Exits with error messages if invalid arguments are given
if re.match(r'\d(?:\.\d{1,4}){3}\/FF\/\d{3,6}', args.funfam) is None:
	sys.exit('Error: -f/--funfam argument must be in the form "C.A.T.H/FF/[FunFam ID]" e.g. 3.40.640.10/FF/63207')
if args.type not in ['strong', 'weak', 'all']:
	sys.exit('Error: -t/--type argument must be "strong", "weak" or "all"')
if args.length not in ['fullgenes', 'genefrags']:
	sys.exit('Error: -l/--length argument must be "fullgenes" or "genefrags"')

# Constructs input file paths based on command-line arguments
	# (Fasta of genomescan run 1 hits is a small subset of all sample genes - therefore much quicker to process)
fasta = {sample: './output/' + sample + '_genes_run1_hits.faa' for sample in args.samples}
crh = {sample: './output/' + args.funfam.split('/FF/')[0] + '/' + sample + '_hits_' + args.type + args.length + '.crh' for sample in args.samples}

# Compiles regular expressions with which to scan the '.crh' file format
p_gene_id = re.compile(r'^(\w+?[+-])')
p_ff = re.compile(r'^\S+?\s(.+?\/FF\/\d+)')
p_dom_bounds = re.compile(r'^(?:\S+?\s){4}(\d+?-\d+?)\s')
p_eval = re.compile(r'^(?:\S+?\s){5}(.+?)\s')


for sample in args.samples:

	# Reads the required data (gene ID, FunFam ID, and domain boundaries) into a nested list (bad practice for heterogenous data, I know)
	with open(crh[sample], 'r') as f:
		crh_data = list(filter(None, f.readlines()))
	crh_list = []
	for hit in crh_data:
		# Only includes hits in the FunFam of interest
		if p_ff.search(hit).group(1) == args.funfam:
			crh_list.append([p_gene_id.search(hit).group(1), p_dom_bounds.search(hit).group(1), \
							p_ff.search(hit).group(1), p_eval.search(hit).group(1)])

	# Reads fasta file into a dict where key=header and value=sequence
	with open(fasta[sample], 'r') as f:
		genes_data = list(filter(None, f.read().splitlines()))
	gene_headers = []
	for line in itertools.islice(genes_data, 0, None, 2):
		gene_headers.append(line.lstrip('>'))
	gene_seqs = []
	for line in itertools.islice(genes_data, 1, None, 2):
		gene_seqs.append(line)
	genes_dict = {}
	for i in range(0, len(gene_headers)):
		genes_dict[gene_headers[i]] = gene_seqs[i]

	# Creates a new 'fasta dict' containing only hits of interest
	# The first letter of the sample name is appended to the start of the gene ID
	# Domain boundaries, FunFam ID, and E-value are also added to the header
	# The domain sequences are extracted from the full sequences
	dom_dict = {}
	for hit in crh_list:
		if hit[0] in gene_headers:  # Filters out sequences that aren't in the crh list
			header = hit[0].replace('k141', '>' + sample[0]) + '_' + hit[1] + ' ' + hit[2] + ' ' + hit[3]
			dom_start = int(hit[1].split('-')[0])
			dom_end = int(hit[1].split('-')[1])
			seq = genes_dict[hit[0]][dom_start-1:dom_end]
			dom_dict[header] = seq

	# If the '-j' command-line option is given, the fasta sequences for all samples are appended to one output file
	if args.joined:
		file_path = './output/' + args.funfam.split('/FF/')[0] + '/' + args.type + args.length + '_' + args.funfam.split('/FF/')[1] + '.faa'
		file_mode = 'a'
	else:
		file_path = './output/' + args.funfam.split('/FF/')[0] + '/' + sample + '_' + args.type + args.length + '_' + args.funfam.split('/FF/')[1] + '.faa'
		file_mode = 'w'
	with open(file_path, file_mode) as out:
		for dom in dom_dict:
			out.write(dom + '\n')
			out.write(dom_dict[dom] + '\n')
