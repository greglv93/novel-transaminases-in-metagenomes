#!/bin/bash

'''
Iterates through a list of HMM files and effectively runs an all vs all profile comparison
Prerequisite: run ../Scripts/extract_funfamHMMs.py to create the necessary HMM library files
Next step: use ./clean_prc_table.py to create a cleaned table ready to be imported into Cytoscape for viewing of the network
'''

# creates output folder in current directory
mkdir ./output

# Runs a one vs all profile comparison for all each of the 157 profile HMMs of interest
# (this may take ~30 minutes minutes since it's calculating 157^2 comparisons)
for file in ../funfam_hmms/hmmer2/*.hmm
do 
	# ${file:22} is file name, not including the directory path (the first 22 characters)
	./prc-1.5.6-linux-x86_64 -algo forward -mode global-global $file ../funfam_hmms/hmmer2/TAm-funfam-hmm2-v4_1_0-list.lib ./output/${file:22}
done

# Concatenates the resulting output files into one large 'all vs all' output
cat ./output/* > ./global_global_PRC
