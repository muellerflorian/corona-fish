
# TODO: include missing blasts (Mycobacterium tuberculosis, Influenza A)
# TODO: important or not to include genetic variability within virus?

# %% Imports
from pathlib import Path
import pandas as pd
from covfish import probe_designer
import importlib


# %% Folders for data and results
path_data = Path(Path.cwd() / '..' / '..' / 'data').resolve()
path_blast = path_data / 'fasta' / 'Probes__cov-2' / 'blast'
path_probes = path_data / 'fasta' / 'Probes__cov-2'
path_results = path_blast / 'results_details'

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


# %% Read probe summary and add other blast results

# >>>  DF with probe list
file_summary = path_probes / 'Probes__cov-2_ALL_summary.txt'
probes_summary = pd.read_csv(file_summary, sep='\t')
probes_summary['qseqid'] = probes_summary['ProbesNames'] + '--' + probes_summary['theStartPos'].astype('str') + '--' + probes_summary['theEndPos'].astype('str')

# >>> Loop over all blast results and join 
probes_summary_joined = probes_summary
blast_all = blast_web_all + blast_local_all 
for blast_add in blast_all:
    probes_summary_joined = pd.merge(left=probes_summary_joined, right=blast_add, how='left', left_on='qseqid', right_on='qseqid')

# >>> Save file to csv
file_save = path_probes / 'Probes__cov-2_ALL__with_blast.csv'
probes_summary_joined.fillna('NAN').to_csv(file_save, sep=',')   
