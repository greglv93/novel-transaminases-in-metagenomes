#!/usr/bin/python3

'''
Some data wrangling that results in a table printed to the terminal window.
This can be copied into a word processor and converted into a table.
'''


import pandas as pd


pd.options.display.float_format = '{:,.0f}'.format
pd.set_option('max_colwidth', 80)
pd.set_option('max_rows', 150)
pd.set_option('display.width', 1000)

pil_wk = pd.read_csv('./output/3.40.640.10/pilluana_stats_weakall', sep='\t')
pil_wk.columns = ['SF','FF','pilluana']
mar_wk = pd.read_csv('./output/3.40.640.10/maras3_stats_weakall', sep='\t')
mar_wk.columns = ['SF','FF','maras3']
all_wk = pd.merge(pil_wk, mar_wk, how='outer', on=['SF','FF'])
all_wk = all_wk.drop('SF', axis=1)

ff_table = pd.read_csv('./3.40.640.10-v4_2_0-funfams-&EC=2.6.1.-_swissprot', sep='|')
ff_table = ff_table.drop(['size','dops','structure','%_2.6.1.-','%_EC4_mode'], axis=1)
ff_table = ff_table.drop(0, axis=0)
ff_table = ff_table.rename(columns={'id':'FF', 'funfam_name':'Name'})

for i in range(0, len(all_wk)):
	all_wk.loc[i,'FF'] = str(all_wk.loc[i,'FF'])
	if isnan(all_wk.loc[i,'pilluana']):
		all_wk.loc[i,'pilluana'] = 0
	else:
		all_wk.loc[i,'pilluana'] = int(all_wk.loc[i,'pilluana']) # still remains a float for some reason
	if isnan(all_wk.loc[i,'maras3']):
		all_wk.loc[i,'maras3'] = 0
	else:
		all_wk.loc[i,'maras3'] = int(all_wk.loc[i,'maras3'])

all_wk['FF'] = all_wk.FF.astype(object)  # doesn't work for row 0 for some reason
all_wk.loc[0,'FF'] = '63148'
all_wk = pd.merge(ff_table,all_wk,how='inner',on='FF')

pil_wk_hits = [x for x in all_wk['pilluana']]
mar_wk_hits = [x for x in all_wk['maras3']]
pil_minus_mar = pd.Series([i - j for i,j in zip(pil_wk_hits, mar_wk_hits)])
all_wk['pil-mar'] = pil_minus_mar.values
all_wk = all_wk.sort_values(by='pil-mar', ascending=False)
all_wk = all_wk[['FF','Name','EC4_mode','pilluana','maras3','pil-mar']]
all_wk = all_wk.rename(columns={'FF':'Closest FF'})

print(all_wk)
print(all_wk.to_csv(sep='|'))
