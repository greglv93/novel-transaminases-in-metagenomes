#!/usr/bin/python3

'''
Creates a FunFam table for a given CATH Superfamily.
Scrapes information from the CATH website (alignments page and stockholm files).
The stockholm files are used to calculate the modal EC annotations for each FunFam, 
	and what percentage of reviewed (SwissProt) sequences are from the desired enzyme class.
The output table is created with only a subset of the FunFams, based on user-defined parameters
	e.g. dominant EC == 2.6.1, DOPS>=70, contains >1 structure, etc
This table can then be merged with downstream processes (such as the metagenome hits in each FunFam)
Prerequisite: SWISSPROT_IDS_PATH, downloaded from UniProtKB
'''

import re
import os
import urllib.request
import itertools
import sys
import argparse

SWISSPROT_IDS_PATH = '/SwissProt_2.6.1.-_26Jun_ids'

# Creates parser object and defines optional arguments
parser = argparse.ArgumentParser(description='Writes a FunFam table for a given CATH Superfamily', add_help=True)
parser.add_argument('-s', '--superfam', dest='superfam', required=True, help='CATH ID e.g. 3.40.640.10')
parser.add_argument('-v', '--version', dest='version', default='latest', help='enter CATH version in format "v4_1_0"')
parser.add_argument('-e', '--ec', dest='ec', required=True, help='enter desired 3-level EC number e.g. 2.6.1')
parser.add_argument('-r', '--reviewed', action='store_true', help='if -r option is given, only reviewed (SwissProt) sequences are used \
	for EC calculations, else, all FunFam members are used')
parser.add_argument('-f', '--function', action='store_true', help='if -f option is given, only FunFams with >=50%% of the desired EC \
	are used for the output table')
parser.add_argument('-3', '--structure', action='store_true', help='if -3 (3d) option is given, only FunFams with at least one 3d structure are included')
parser.add_argument('-d', '--dops', action='store_true', help='if -d option is given, only FunFams with DOPS>70 (high diversity) are included')
# Reads in command-line arguments and adds the attributes to the parser object
args = parser.parse_args()

# Exits with error messages if invalid arguments are given
if re.match(r'\d(\.(\d){1,4}){3}', args.superfam) is None:
	sys.exit('Error: a correct CATH code was not given in the first argument')
if re.match(r'\d\.\d{1,2}\.\d{1,3}', args.ec) is None:
	sys.exit('Error: -e/--ec argument must be a valid 3-level EC code')


# Scrapes the funfam table info from the CATH webpage
url = 'http://www.cathdb.info/version/' + args.version + '/superfam/' + args.superfam + '/alignments'
with urllib.request.urlopen(url) as response:
   html = str(response.read())
# Splits html text so that each list item starts with a FunFam number
html = html.split('data-funfam-number=')[1:]

# Regular expressions to find all the FunFam data from the html string
p_ffid = re.compile(r'^"(\d+?)"')
p_name = re.compile(r'data-funfam-name="([\s\S]+?)"')
p_dops = re.compile(r'data-funfam-dops="([\d\.]+?)"')
p_size = re.compile(r'data-funfam-members="([\s\S]+?)"')
p_3d = re.compile(r'\<div class="label label-3D"\>3D\<\/div\>')
p_rep = re.compile(r'data-funfam-rep="([\s\S]+?)"')
# Extracts the FunFam data into a dict, including filtering based on optional arguments
funfam_dict = {}
for funfam in html:
	ffid = p_ffid.search(funfam).group(1)
	name = p_name.search(funfam).group(1)
	dops = p_dops.search(funfam)
	if dops is not None:
		dops = float(dops.group(1))
	else:
		continue  # DOPS = 0
	if args.dops:
		if dops < 70:
			continue
	size = int(p_size.search(funfam).group(1))
	structure = p_3d.search(funfam)
	if structure is not None:
		structure = p_rep.search(funfam).group(1)
	else:
		if args.structure:
			continue
		else:
			structure = '-'
	funfam_dict[ffid] = [name, size, dops, structure]


# If only SwissProt sequences are to be used in the EC calculations, a file of all SwissProt IDs is loaded into memory
if args.reviewed:
	SwissProt_ID_path = os.getcwd() + SWISSPROT_IDS_PATH
	with open(SwissProt_ID_path, 'r') as f:
		SwissProt_IDs = f.read().splitlines()


# Patterns to find accession code and EC annotation in stockholm files
p_ac = re.compile(r'^(\S+?)\\n')
p_ec4 = re.compile(r'DR\sEC;((\s\d+?\.\d+?\.\d+?\.\d+?;)+)')
# Gets stockholm file for each FunFam and makes 3 calculations: % target EC, modal EC4, and % modal EC4,
	# which are appended to the funfam dict, unless the -f argument is given and target EC < 50%, in which case the FunFam is removed
for funfam in funfam_dict:
	sto_url = 'http://www.cathdb.info/version/' + args.version + '/superfam/' + superfam + '/funfam/' + funfam + '/files/stockholm'
	with urllib.request.urlopen(sto_url) as response:
		sto = str(response.read())
	members = sto.split(' AC ')[2:]
	target_ec3_count = 0
	total_ec_count = 0
	ec4_list = []
	for member in members:
		seq_ac = p_ac.search(member)
		if seq_ac is None:
			continue
		if args.reviewed:
			seq_ac = seq_ac.group(1)
			if seq_ac not in SwissProt_IDs:
				continue
		ec4 = p_ec4.search(member)
		if ec4 is not None:
			total_ec_count += 1
			ec4 = ec4.group(1)
		else:
			continue
		if args.ec in ec4:
			target_ec3_count += 1
		ec4 = ec4.split(';')[:-1]
		for code in ec4:
			code = code.lstrip(' ')
			ec4_list.append(code)
	if target_ec3_count == 0:
		target_ec3_perc = 0
	else:
		target_ec3_perc = round((target_ec3_count/total_ec_count)*100, 2)
	if args.function and target_ec3_perc < 50:
		funfam_dict.pop(funfam, None)
		continue
	else:
		funfam_dict[funfam].append(str(target_ec3_perc))
	if len(ec4_list) == 0:
		ec4_mode = '-'
		ec4_mode_perc = 0
	else:
		ec4_mode = max(set(ec4_list), key=ec4_list.count)
		ec4_mode_perc = round((ec4_list.count(ec4_mode)/total_ec_count)*100, 2)
	funfam_dict[funfam].append(ec4_mode)
	funfam_dict[funfam].append(str(ec4_mode_perc))


# Configures some file name variables based on the optional arguments given
if args.dops:
	dops_file_name = 'dops=70'
else:
	dops_file_name = ''
if args.structure:
	struc_file_name = '&3d=yes'
else:
	struc_file_name = ''
if args.function:
	func_file_name = '&EC=' + args.ec + '.-'
else: func_file_name = ''
if args.reviewed:
	swissprot_or_uniprot = 'swissprot'
else:
	swissprot_or_uniprot = 'uniprot'

# Writes a flat text table with separator='|'
out_path = cwd + '/' + superfam + '-' + args.version + '-funfams-' + dops_file_name + struc_file_name + func_file_name + '_' + swissprot_or_uniprot
with open(out_path, 'w') as out:
	out.write('id|size|dops|structure|%_' + args.ec + '.-|EC4_mode|%_EC4_mode|funfam_name' + '\n' + '-'*70 + '\n')
	for funfam in funfam_dict:
		if args.ec:
			if float(funfam_dict[funfam][4]) < 50:
				continue
		out.write(funfam + '|' + str(funfam_dict[funfam][1]) + '|' + str(funfam_dict[funfam][2]) + '|')
		out.write(funfam_dict[funfam][3] + '|' + funfam_dict[funfam][4] + '|')
		out.write(funfam_dict[funfam][5] + '|' + funfam_dict[funfam][6] +  '|' + funfam_dict[funfam][0] + '\n')
