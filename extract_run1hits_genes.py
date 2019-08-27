#!/usr/bin/python3

'''
Takes in .crh file from genomescan run 1 and extracts fasta sequences of all genes that contain hits.
'''

import argparse
import re
import itertools
import os


# Creates parser object and defines optional arguments
parser = argparse.ArgumentParser(description='Takes in .crh file from genomescan run 1 and extracts fasta sequences of all genes that contain hits', \
	add_help=True)
parser.add_argument('-i', '--input', dest='genes_path', help='Enter path for all genes fasta (the FragGeneScan output)')
parser.add_argument('-c', '--crh', dest='crh_path', help='Enter path for .crh file from genomescan run 1')
parser.add_argument('-o', '--out', dest='out_path', help='Enter output fasta path')
# Reads in command-line argument and adds the supfam attribute to the parser object
args = parser.parse_args()

# Reads in fasta sequence file containing all genes
with open(args.genes_path, 'r') as f:
	genes_data = f.read().splitlines()
genes_data = list(filter(None, genes_data))
gene_headers = []
for line in itertools.islice(genes_data, 0, None, 2):
	gene_headers.append(line)
gene_seqs = []
for line in itertools.islice(genes_data, 1, None, 2):
	gene_seqs.append(line)
genes_dict = {}
for i in range(0, len(gene_headers)):
	genes_dict[gene_headers[i]] = gene_seqs[i]

# Reads in crh output from genomescan run 1
with open(args.crh_path, 'r') as f:
	crh = f.read().splitlines()
crh = list(filter(None, crh[2:]))
hits_genes = []
p_id = re.compile(r'^(\w+?[+-])')
for hit in crh:
	gene_id = p_id.search(hit).group(1)
	gene_id = '>' + gene_id
	if gene_id in hits_genes:
		continue
	else:
		hits_genes.append(gene_id)

# Writes output fasta using only gene IDs with hits in genomescan run 1
with open(args.out_path, 'w') as out:
	for x in hits_genes:
		out.write(x + '\n' + genes_dict[x] + '\n')

# gets rid of last '\n'
with open(out_path, 'rb+') as out:
	out.seek(-1, os.SEEK_END)
	out.truncate()
