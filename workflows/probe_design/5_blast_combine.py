
# TODO: include missing blasts (Mycobacterium tuberculosis, Influenza A)
# TODO: important or not to include genetic variability within virus?

# %% Imports
import sys
from pathlib import Path
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# %% Read different blast results

path_blast = Path(Path.cwd() / '..' / '..' / 'data' / 'fasta' / 'Probes__cov-2' / 'blast').resolve()
path_probes = Path(Path.cwd() / '..' / '..' / 'data' / 'fasta' / 'Probes__cov-2').resolve()

# Local blast against corona (separator is \t)
file_blast = path_blast / 'corona_family.txt'
blastn_corona_family = pd.read_csv(file_blast, sep='\t', header=None,
                           names = ['qseqid','cf_sseqid','cf_pident','cf_length','cf_mismatch','cf_gapopen','cf_qstart','cf_qend','cf_sstart','cf_send','cf_evalue','cf_bitscore'])

# Internet blast against human  (separator is ')
file_blast =  path_blast / 'blastn_human_transcriptome.csv'
blastn_human = pd.read_csv(file_blast, sep=',', header=None,
                           names = ['qseqid','h_sseqid','h_pident','h_length','h_mismatch','h_gapopen','h_qstart','h_qend','h_sstart','h_send','h_evalue','h_bitscore'])

# Remove NC entries from human blast
# TODO: MAYBE KEEP ONLY TRANSCRIPT BASED STUFF: NM	mRNA,  NR ncRNA,  XM predicted mRNA model, XR predicted ncRNA model
refseq_strings = ['NM_','XM_','NR','XR_']
blastn_human = blastn_human[blastn_human.h_sseqid.str.contains('|'.join(refseq_strings ))]

# Get the entries with the highest bitscore
df_sorted_1 = blastn_corona_family.sort_values(by='cf_bitscore')
corona_best_hit = df_sorted_1.drop_duplicates('qseqid', keep='last').sort_index(0)

df_sorted_2 = blastn_human.sort_values(by='h_bitscore')
human_best_hit = df_sorted_2.drop_duplicates('qseqid', keep='last').sort_index(0)

corona_best_hit.drop(columns= ['cf_pident','cf_gapopen','cf_qstart','cf_qend','cf_sstart','cf_send','cf_evalue'], inplace=True)
human_best_hit.drop(columns= ['h_pident','h_gapopen','h_qstart','h_qend','h_sstart','h_send','h_evalue'], inplace=True)


# >>>> DF with list of probes
file_summary = path_probes / 'Probes__cov-2_ALL.txt'
probes_summary = pd.read_csv(file_summary, sep='\t')

# Create matching column name with blast results
probes_summary['qseqid'] = probes_summary['ProbesNames'] + '--' + probes_summary['theStartPos'].astype('str') + '--' + probes_summary['theEndPos'].astype('str')

probes_summary = pd.merge(left=probes_summary, right=corona_best_hit, how='left', left_on='qseqid', right_on='qseqid')
probes_summary = pd.merge(left=probes_summary, right=human_best_hit, how='left', left_on='qseqid', right_on='qseqid')


# %% Save as csv
file_save = path_probes / 'Probes__cov-2_ALL_with_blast.csv'
probes_summary.fillna('NAN').to_csv(file_save, sep=',')


# %% Show some results
# Legends thanks to this: https://stackoverflow.com/questions/29096632/getting-legend-in-seaborn-jointplot/29909033

# Blast against coronavirus
file_save = path_blast / 'blast_cf_mismatch.png'
hexplot = sns.jointplot("cf_length", "cf_mismatch", data=probes_summary, kind="hex")       
plt.subplots_adjust(left=0.2, right=0.8, top=0.8, bottom=0.2)  # shrink fig so cbar is visible
cbar_ax = hexplot.fig.add_axes([.85, .25, .05, .4])  # x, y, width, height
plt.colorbar(cax=cbar_ax)
plt.savefig(file_save,dpi =300)

# Blast against human transcriptome
file_save = path_blast / 'blast_h_mismatch.png'
hexplot = sns.jointplot("h_length", "h_mismatch", data=probes_summary, kind="hex")       
plt.subplots_adjust(left=0.2, right=0.8, top=0.8, bottom=0.2)  # shrink fig so cbar is visible
cbar_ax = hexplot.fig.add_axes([.85, .25, .05, .4])  # x, y, width, height
plt.colorbar(cax=cbar_ax)
plt.savefig(file_save,dpi =300)

# %%

