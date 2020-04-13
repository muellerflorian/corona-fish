
# %% Imports
import pandas as pd
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# %% Specify folders and files
path_probes = Path(Path.cwd() / '..' / '..' / 'data' / 'probes' / 'Probes__cov-2').resolve()
file_results = path_probes / 'Probes__cov-2_ALL__with_blast_FLAPY__coverage.csv'
probes_summary_load = pd.read_csv(file_results)
probes_summary_load = probes_summary_load.fillna(0)

# %% GC and PNAS filters
query_probes = 'GCFilter==1 & NbOfPNAS>2'
n_probes = probes_summary_load.query(query_probes).shape[0]
print(f'Query [{query_probes}] yields {n_probes} probes')

query_cov_align = 'cov2_align_perc > 98'
n_probes = probes_summary_load.query(query_cov_align).shape[0]
print(f'Query [{query_cov_align}] yields {n_probes} probes')

query_both = query_probes + '&' + query_cov_align
n_probes = probes_summary_load.query(query_both).shape[0]
print(f'Query [{query_both}] yields {n_probes} probes')

# %%  Get all alignment length in file
blast_align_max = 22

col_names = list(probes_summary_load.columns.values)
cols_length = [col_name for col_name in col_names if ('_length' in col_name) and (not 'cov-2' in col_name)  ]

for col_length in cols_length:
    col_mismatch = col_length.replace('length','mismatch')
    query = f'{col_length}-{col_mismatch}<={blast_align_max}'
    n_probes = probes_summary_load.query(query).shape[0]
    print(f'Query [{query}] yields {n_probes} probes')

# %% Query all
query_align_length = ''.join(map(lambda x: x + '-' + x.replace('length', 'mismatch') + f'<={blast_align_max} & ', cols_length)) 
query_align_length = query_align_length[:-2]  #
n_probes = probes_summary_load.query(query_align_length).shape[0]
print(f'Query [COMBINED alignment length] yields {n_probes} probes')

query_all = query_align_length + '&' + query_probes + '&' + query_cov_align
probes_query = probes_summary_load.query(query_all)
n_probes = probes_query.shape[0]
print(f'Query [ALL restrictions] yields {n_probes} probes')


# %%  Select probes with highest coverage
probes_query = probes_query.nlargest(96, columns=['ngs_cov'])

# %%
path_save = path_probes / 'query' / f'blast_align_max_{blast_align_max}'
if not path_save.is_dir():
    path_save.mkdir(parents=True)

# Save query
file_save = path_save / 'query_string.txt'
with open(file_save, "w") as text_file:
    print(f"{query_all}", file=text_file)

# Plot queried probe list
file_save = path_save / 'cov2_probes_query.csv'
probes_query.to_csv(file_save, sep=',', index=False) 

# Save plot with probe positions
plt.figure(figsize=(10,1))
plt.plot(probes_query['theStartPos'], np.ones((probes_query.shape[0],1)), '|', color='black')
plt.tight_layout()
file_save = path_save / 'probes_positions.png'
plt.savefig(file_save, dpi=300)

# Save plot with probe positions
plt.figure(figsize=(10,1))
plt.plot(probes_query['theStartPos'], np.ones((probes_query.shape[0],1)), '|', color='black')
plt.tight_layout()
file_save = path_save / 'probes_positions.png'
plt.savefig(file_save, dpi=300)

# %%
