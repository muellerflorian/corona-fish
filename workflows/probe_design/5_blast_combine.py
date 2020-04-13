
# %% Specify which data-set should be analyzed
name_probes = 'Probes__cov-2'  # Genomic probes


# %% Imports
from pathlib import Path
import pandas as pd
from covfish import probe_designer
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import importlib


# %% Folders for data and results
path_data = Path(Path.cwd() / '..' / '..' / 'data').resolve()
path_blast = path_data / 'probes' / f'{name_probes}' / 'blast'
path_probes = path_data / 'probes' / f'{name_probes}'
path_results = path_blast / 'results_details'

file_summary = path_probes / f'{name_probes}_ALL_summary.txt'

if not path_results.is_dir():
    path_results.mkdir()

# %%  Read blast results
importlib.reload(probe_designer)

# >>> Import web blast
blast_web_all = probe_designer.blast_summarize(path_blast=path_blast / 'web',
                                               file_identifier=path_data / 'blast_identifiers.json',
                                               sep=',', 
                                               ext='csv',
                                               path_results=path_results,
                                               refseq_keep=['NM_', 'XM_', 'NR', 'XR_'])


# >>> Import lcoal blast against other viruses
importlib.reload(probe_designer)
blast_local_all = probe_designer.blast_summarize(path_blast=path_blast / 'local',
                                                 file_identifier=path_data / 'blast_identifiers.json',
                                                 sep='\t',
                                                 ext='txt',
                                                 path_results=path_results,
                                                 refseq_keep=None)

#  >>> import blast against cov-2 alignment

# Count how many genoems (identified by >)
file_cov = path_data / 'genomes' / 'cov2' / 'cov2_aligned.fasta'
genomes_cov2 = file_cov.read_text()
n_genomes = genomes_cov2.count('>')

# Open blast results
names_col = ['qseqid', 'sseqid', 'qlen', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'hsend', 'evalue', 'bitscore']
file_blast = path_blast / 'local_cov2' / 'cov2_aligned.txt'
blast_results = pd.read_csv(file_blast, sep='\t', header=None, names=names_col)
blast_results['align_diff'] = blast_results['qlen'] - blast_results['length'] + blast_results['mismatch'] 

# Simplify df
blast_results = blast_results.filter(items=['qseqid', 'sseqid', 'align_diff'])
# Remove potential duplicates and only heep the ones with highest alignment difference
blast_results.sort_values('align_diff', ascending=False).drop_duplicates(['qseqid','sseqid'], inplace=True)
# Keep only blast hits with 0 or 1 mismatch over the entire probe sequence
blast_covid_alignment = blast_results.query('align_diff <=1').groupby('qseqid').count()
blast_covid_alignment['align_diff'] = blast_covid_alignment['align_diff'] * 100.0 / n_genomes
blast_covid_alignment.rename(columns={'sseqid': f'cov2_align_N',
                                      'align_diff': f'cov2_align_perc'}
                             ,inplace=True)

# >>>>  Read probe summary and add other blast results

# >>>  DF with probe list
probes_summary = pd.read_csv(file_summary, sep='\t')
probes_summary['qseqid'] = probes_summary['ProbesNames'] + '--' + probes_summary['theStartPos'].astype('str') + '--' + probes_summary['theEndPos'].astype('str')

# >>> Loop over all blast results and join 
probes_summary_joined = probes_summary
blast_all = [blast_covid_alignment] + blast_web_all + blast_local_all 
for blast_add in blast_all:
    probes_summary_joined = pd.merge(left=probes_summary_joined, right=blast_add, how='left', left_on='qseqid', right_on='qseqid')

# >>> Save file to csv
file_save = path_probes / f'{name_probes}_ALL__with_blast_FLAPY.csv'
probes_summary_joined.to_csv(file_save, sep=',', index=False)
print(f'\n\nSUMMARY saved as {file_save}')

# %% Plot results of alignment against cov-2

# >>> Number of probe squences that match well all consensus sequences
plt.figure(figsize=(5, 4))
sns.distplot(blast_covid_alignment[f'cov2_align_perc'], 
             bins = np.arange(90,100.05,0.25),
             kde=False, rug=False)
plt.title(f'Percentage of alignment\nwith {n_genomes} genomes for each oligo (n={blast_covid_alignment.shape[0]})')
plt.tight_layout()
file_save = path_results / 'cov2-alignment_hist.png'
plt.savefig(file_save, dpi=300)

# >> Cummulative histogram
plt.figure(figsize=(5, 4))
kwargs = {'cumulative': True}
sns.distplot(blast_covid_alignment[f'cov2_align_perc'], hist_kws=kwargs, kde_kws=kwargs)
plt.title(f'Percentage of alignment\nwith {n_genomes} genomes for each oligo (n={blast_covid_alignment.shape[0]})')
file_save = path_results / 'cov2-alignment_cdf.png'
plt.savefig(file_save, dpi=300)
# %%
