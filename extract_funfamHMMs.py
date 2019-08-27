#!/usr/bin/python3

'''
Takes in a large HMMER3 profile library (all FunFams in CATH v4.1) and extracts profiles of a desired subset of CATH Superfamilies.
Writes one HMMER3 library file with that subset, and also separate '.hmm' files for each individual profile.
Also converts the profiles of interest to HMMER2 format and writes the same output as for the HMMER3 files.
(Some programs, such as Profile Comparer, only take HMMER2 files as input)
Prerequisite: there should be a folder ../funfam_hmms/ (relative to this script)
Next step: use ../PRC/iterative_PRC.sh to run a profile comparison between all HMMs extracted here
'''

import argparse
import re
import os
import subprocess


# Creates parser object and defines optional arguments
parser = argparse.ArgumentParser(description='Takes in HMMER3 profile library for all FunFams in CATH v4.1 \
									and creates a smaller profile library of profiles in a given CATH superfamily', add_help=True)
parser.add_argument('-s', '--supfams', dest='supfams', required=True, nargs='+', help='list of one or more superfamilies e.g. 3.40.640.10')
# Reads in command-line arguments and adds the 'supfams' attribute to the parser object
args = parser.parse_args()
# Exits with error messages if invalid arguments are given
for supfam in args.supfams:
	if re.match(r'\d(?:\.\d{1,4}){3}', supfam) is None:
		sys.exit('Error: -s/--supfams argument must be in the form "C.A.T.H" e.g. 3.40.640.10')


# Pattern that matches the FunFam code in the HMMER3 file format
# group(1) is the full code, group(2) matches only the Superfamily code
p_id = re.compile(r'NAME\s+((\S+?)\/FF\/\d+)')


hmm = ''
tam_hmms = []
# Reads in the full HMMER3 file with profile HMMs of all 92,882 HMMs
with open('../genomescan/data/funfam-hmm3-v4_1_0.lib', 'r') as infile:
	# Iterates through each line and only adds the HMM if in the subset of interest
	for line in infile:
		# Adds each line to the 'current' HMM
		if '//' not in line:
			hmm += line
		# When end of current HMM is reached, it is added to the list if in the subset of interest
		elif '//' in line:
			cathid = p_id.search(hmm).group(2)
			for supfam in args.supfams:
				if supfam == cathid:
					tam_hmms.append(hmm)
			# Resets the current HMM variable
			hmm = ''

# Creates a new HMMER3 library file with only the subset of profiles of interest
# This will be used in the first run of the genomescan program
with open('../genomescan/data/TAm-funfam-hmm3-v4_1_0.lib', 'w') as out:
	out.write('//\n'.join(tam_hmms))
	out.write('//')

# Writes a separate file for each HMM - these are needed as input for the PRC program
for hmm in tam_hmms:
	hmm_name = p_id.search(hmm).group(1).replace('/', '.')
	with open('../funfam_hmms/' + hmm_name + '.hmm', 'w') as out:
		out.write(hmm)
		out.write('//')

# Converts individual HMM files to HMMER2 format afterwards as PRC doesn't accept HMMER3 files
p_hmmfile = re.compile(r'.+\.hmm$')  # used to check if a file ends in '.hmm'
hmms_path = os.getcwd().rstrip('Scripts') + 'funfam_hmms/'
files = os.listdir(hmms_path)
subprocess.check_call(['mkdir', '../funfam_hmms/hmmer2'])
for file in files:
	if p_hmmfile.search(file) is not None:
		with open(hmms_path + 'hmmer2/' + file, 'w') as f:
			subprocess.check_call(['../genomescan/bin/hmmer3/hmmconvert', '-2', hmms_path + file], stdout=f)

# Writes a different library file: a simple text file listing the file paths of the individual HMM profiles
# This is also needed as input for the PRC program
with open('../funfam_hmms/hmmer2/TAm-funfam-hmm2-v4_1_0-list.lib', 'w') as out:
	files = os.listdir(hmms_path + 'hmmer2/')
	for file in files:
		if p_hmmfile.search(file) is not None:
			out.write(hmms_path + 'hmmer2/' + file)
			out.write('\n')
