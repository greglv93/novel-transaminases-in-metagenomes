#!/usr/bin/python3

import pandas as pd
import math
import seaborn as sn
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors


# reads FunFam IDs into a DataFrame
df = pd.read_csv('./inc_thresh/3.40.640.10', sep='\t', header=None)  # just a file that happens to have a list of all the FFs in the SF in one column
df = df.rename(columns={0:'FF', 1:'drop1', 2:'drop2'})
PF_FF_df = df.drop(['drop1', 'drop2'], axis=1)  # only the 'FF' column remains - this will be set as the index after merging with the mapping data

# Pfam IDs and corresponding mapping files
pfam_fams = ['PF00155', 'PF00202', 'PF00266', 'PF01041', 'PF12897']  # doesn't include PF01063
pfam_ff_mappings = {pfam: '../mappings/' + pfam + '_crh_catalytic_funfam_composition' for pfam in pfam_fams}  # dict with file names & locations

# adds each pfam mapping to the funfam dataframe
for pfam in pfam_ff_mappings:
	df = pd.read_csv(pfam_ff_mappings[pfam], sep='\t', header=None)
	df.columns = [pfam, 'FF']
	df = df[['FF', pfam]]
	for i in range(0, len(df)):
		df.loc[i, 'FF'] = int(df.loc[i, 'FF'].split('FF/')[-1])  # changed to int to make dtype same as merging df
	PF_FF_df = pd.merge(PF_FF_df, df, how='outer', on='FF')

# make 'FF' the index
PF_FF_df = PF_FF_df.set_index('FF')
# make sure the columns are in the desired order (as specified above)
PF_FF_df = PF_FF_df[pfam_fams]

# algorithm that sorts rows column by column, limited to values that are highest in that column
	# this should give some order to the visual presentation of the heatmap
PF_FF_df_copy = PF_FF_df.copy()  # copies original DataFrame, will be used as a reference and emptied as rows are copied
PF_FF_df = pd.DataFrame(columns=pfam_fams)  # empties the original DataFrame; it will be filled with rows from the copy
for pfam in PF_FF_df_copy:
	PF_FF_df_copy = PF_FF_df_copy.sort_values(by=pfam, ascending=False)  # sorts rows by values in current column
	for FF in list(PF_FF_df_copy.index.values):  # iterates through rows
		if (PF_FF_df_copy.loc[FF,pfam] > 0) and (PF_FF_df_copy.loc[FF,pfam] == PF_FF_df_copy.loc[FF,].max()):
			PF_FF_df = PF_FF_df.append(PF_FF_df_copy.loc[FF])  # appends the row
			PF_FF_df_copy = PF_FF_df_copy.drop(FF, axis=0)  # drops row from reference DataFrame

# adds the remaining FFs with no hits in any Pfam
PF_FF_df = pd.concat([PF_FF_df, PF_FF_df_copy])

# writes the DataFrame to an output file	
PF_FF_df.to_csv('./output/PF_FF_mapping', sep='\t')

# preparing the DataFrame for heatmap
PF_FF_df = PF_FF_df.applymap(math.log10)  # converts values to log10
PF_FF_df = PF_FF_df.fillna(0)  # replaces all NaN with 0

# creates figure and axes objects
fig, ax = plt.subplots(figsize=(20, 4))

# plots heatmap with transposed dataframe
ax = sn.heatmap(
	PF_FF_df.T,
	cbar_kws={'ticks': [1,2,3,4],
	'label': 'log10(sequences)'})

# rotates y tick labels
for item in ax.get_yticklabels():
	item.set_rotation(0)
# changes font size
matplotlib.rcParams.update({'font.size': 12})
# adds axes lables
plt.xlabel('CATH FunFams in Superfamily 3.40.640.10')
plt.ylabel('Pfam families')
# removes plot ticks and x-axis labels (too crowded)
plt.tick_params(
	axis='x',          # changes apply to the x-axis
	which='both',      # both major and minor ticks are affected
	bottom=False,      # ticks along the bottom edge are off
	top=False,         # ticks along the top edge are off
	labelbottom=False) # labels along the bottom edge are off

plt.savefig('./PF_FF_heatmap_log10.png', dpi=100)
plt.show()


# TO DO:
	# slight inaccuracy - log10(1) = 0 and NaN is converted to 0 (after log conversion), so the heatmap doesn't differentiate between 1 and 0 hits
		# log(x+1) transformation instead?
		# doesn't matter too much since the sorting algorithm works on the original values, so we know which squares in the heatmap are 1
			# for all pfam families except the last, where the 1 values are directly next to the 0 values at the end
	# change legend tick labels from [1,2,3,4] to [10, 100, 1000, 10000] 
		# see 2nd answer on https://stackoverflow.com/questions/11244514/modify-tick-label-text?noredirect=1&lq=1
	# make plot bigger and show all the x-axis (FF) labels - atm only showing ~1/4 
		# or just label/title the axis 'FunFams' without any actual labels
			# just provide the full table in appendices/supplementary work if needed
			# also get rid of x axis ticks
	# change y-axis labels to the format 'Aminotransferase Class-III (PF00202)' and add a y-axis title (Pfam family)
		# or just write in the figure description
	# change plot dimensions: 3-4x larger width relative to height
	# change colour scale legend to show the actual number (reverse the log function) and add more small ticks
	# change background colour to white? (the colour scale would also have to change e.g. darker = higher value)
	# figure description: heatmap gives a good overview of the sequence diversity of each of the Pfams
