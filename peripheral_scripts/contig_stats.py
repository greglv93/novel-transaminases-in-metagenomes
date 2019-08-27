#!/usr/bin/python3

import re

PILLUANA_CONTIGS = "./pilluana.contigs.fa"
MARAS3_CONTIGS = "./maras3.contigs.fa"

p_len = re.compile(r'len=(\d+)')


# read in contig files and split into list
with open(PILLUANA_CONTIGS, 'r') as f:
    pil_contigs = f.read().split('>')

with open(MARAS3_CONTIGS, 'r') as f:
    mar_contigs = f.read().split('>')


# extract length of each contig and put into list of lengths
pil_lengths = []
for pil in pil_contigs:
    if p_len.search(pil) is not None:
        pil_lengths.append(int(p_len.search(pil).group(1)))

mar_lengths = []
for mar in mar_contigs:
    if p_len.search(mar) is not None:
        mar_lengths.append(int(p_len.search(mar).group(1)))


# write lengths into a simple comma-separated file
with open('./pilluana.contigs_lengths', 'w') as out:
    for item in pil_lengths[:-1]:
        out.write(str(item) + ',')
    out.write(str(pil_lengths[-1]))

with open('./maras3.contigs_lengths', 'w') as out:
    for item in mar_lengths[:-1]:
        out.write(str(item) + ',')
    out.write(str(mar_lengths[-1]))


# print out some basic stats
pil_plus_5k_count = 0
for l in pil_lengths:
	if l >= 5000:
		pil_plus_5k_count +=1

print("no. of pilluana contigs over 5kb: %s" % str(pil_plus_5k_count))

pil_perc_plus_5k = float((pil_plus_5k_count*100))/float(len(pil_contigs))

print("percentage of pilluana contigs over 5kb: %.3f" % pil_perc_plus_5k)


mar_plus_5k_count = 0
for l in mar_lengths:
	if l >= 5000:
		mar_plus_5k_count +=1

print("no. of maras3 contigs over 5kb: %s" % str(mar_plus_5k_count))

mar_perc_plus_5k = float((mar_plus_5k_count*100))/float(len(mar_contigs))

print("percentage of maras3 contigs over 5kb: %.3f" % mar_perc_plus_5k)

