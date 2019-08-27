#!/usr/bin/python3

'''
Takes in large fasta file of predicted ORFs (FragGeneScan output) and contig sequences.
Writes plain text file with list of fasta headers of the subset of predicted ORFs that are not truncated by contig boundaries.
Conditions: sequence doesn't start at the first position or end at the last position of the contig
This is quicker than creating a function that needs to be used iteratively to make this assessment for each downstream script
(the large fasta file of all genes would have to be opened and checked every time),
'''

import re
import itertools

# change this for different samples
SAMPLES = ['pilluana', 'maras3']

# Patterns to extract key info from fasta headers
p_id = re.compile(r'^\>(k141_\d+)')
p_len = re.compile(r'\slen=(\d+?)\n')
p_gene_start = re.compile(r'^\>k141_\d+?_(\d+?)_')
p_gene_end = re.compile(r'^\>k141_\d+?_\d+?_(\d+?)_')
p_strand = re.compile(r'^>k141_(?:\d+?_){3}([-\+])')


for sample in SAMPLES:

	#(1)# Gets contig headers and puts id, len as key, value in a dict
	# This will be used in step (3) for comparing the end of the gene to the end of the contig (if the same, the gene has been cut short)
	print('reading in %s contig lengths...' % sample)
	contig_headers = []
	contig_path = '../MG_data/' + sample + '.contigs.fa'
	with open(contig_path, 'r') as f:
		for line in itertools.islice(f, 0, None, 2):
			contig_headers.append(line)
	contig_lengths = {}
	for contig in contig_headers:
		contig_id = p_id.search(contig).group(1)
		contig_len = p_len.search(contig).group(1)
		contig_lengths[contig_id] = contig_len
	print('finished reading in %d contig lengths' % len(contig_lengths))

	#(2)# Gets FGS output (genes) and puts gene start, stop and seq in a dict
	print('reading in FragGeneScan predicted genes for %s...' % sample)
	fgs_path = '../FragGeneScan/results/' + sample + '_genes.faa'
	with open(fgs_path, 'r') as f:
		fgs_data = f.read().splitlines()
	gene_headers = []
	for line in itertools.islice(fgs_data, 0, None, 2):
		gene_headers.append(line)
	gene_dict = {}  # value = (contig_id, start, stop, strand)
	for i in range(0, len(gene_headers)):
		gene_id = gene_headers[i]
		contig_id = p_id.search(gene_headers[i]).group(1)
		gene_start = p_gene_start.search(gene_headers[i]).group(1)  # gene_start is actually the end for '-' strand genes
		gene_end = p_gene_end.search(gene_headers[i]).group(1)  # gene_end is actually the start for '-' strand genes
		gene_strand = p_strand.search(gene_headers[i]).group(1)
		gene_dict[gene_id] = (contig_id, gene_start, gene_end, gene_strand)
	print('finished reading in data for %d genes' % len(gene_dict))

	#(3)# Writes a file with full genes 
	print('writing output file "%s_full_gene_ids" (list of IDs of full-length genes)...' % sample)
	full_length_count = 0
	out_path = './' + sample + '_full_gene_ids'
	with open(out_path, 'w') as out:
		for gene in gene_dict:
			if gene_dict[gene][3] == '+':
				if int(gene_dict[gene][1]) == 1:
					# start of '+' strand geneis likely to be missing
					continue  # skips to next gene
				if int(gene_dict[gene][2]) >= (int(contig_lengths[gene_dict[gene][0]]) - 2):
					# end of '+' strand is missing
					continue  # skips to next gene
				out.write(gene + '\n')
				full_length_count += 1
			elif gene_dict[gene][3] == '-':
				if int(gene_dict[gene][1]) <= 3:
					# end of '-' strand gene is missing
					continue  # skips to next gene
				if int(gene_dict[gene][2]) == int(contig_lengths[gene_dict[gene][0]]):
					# start of '-' strand gene is likely to be missing
					continue  # skips to next gene
				out.write(gene + '\n')
				full_length_count += 1
	print('finished writing file')

	# Prints some information
	total_predicted = len(gene_headers)
	print('total predicted %s genes: %d' % (sample, total_predicted))
	print('number of full-length %s genes: %d (%d %% of total)' 
		% (sample, full_length_count, round( ((full_length_count*100)/total_predicted), 1))
		)

# results:
	# pilluana_genes: 301436 / 1460829 (20%) are full-length
	# maras3_genes: 290016 / 2197977 (13%) are full-length