
# __IMPORTANT__: to run this file, you need to have the coverage data analyzed (resulting in a coverage.csv file)


# %% Specify which data-set should be analyzed
name_probes = 'Probes__cov-2'  # Genomic probes
#name_probes = 'Probes__cov-2--RevComp'  # Replication intermediate

reverse_complement = True
length_genome = 29904

# %% Imports
from pathlib import Path
import pandas as pd

import matplotlib.pyplot as plt

# %% Folders for data and results
path_data = Path(Path.cwd() / '..' / '..' / 'data').resolve()
file_coverage = path_data / 'coverage' / 'coverage.csv'
coverage = pd.read_csv(file_coverage)

path_probes = Path(Path.cwd() / '..' / '..' / 'data' / 'probes' / f'{name_probes}').resolve()
file_results = path_probes / f'{name_probes}_ALL__with_blast_FLAPY.csv'
probes_summary = pd.read_csv(file_results)

path_results = path_probes / 'coverage'
if not path_results.is_dir():
    path_results.mkdir()

# %% Calculate sum of NGS probe coverage for each probe
def calc_coverage(start, end):
    return coverage.query('pos>=@start & pos<=@end')['cov'].median(axis=0)

if not reverse_complement:
    probes_summary['ngs_cov'] = probes_summary.apply(lambda row: calc_coverage(row['theStartPos'], row['theEndPos']), axis=1)
else:
    probes_summary['ngs_cov'] = probes_summary.apply(lambda row: calc_coverage(length_genome-row['theEndPos'], length_genome-row['theStartPos']), axis=1)
    
# Save plot with probe positions
plt.figure(figsize=(10,6))
plt.plot(probes_summary['theStartPos'],probes_summary['ngs_cov'],'o')
plt.xlabel('position')
plt.ylabel('NGS coverage')
plt.tight_layout()
file_save = path_results / 'probes_coverage.png'
plt.savefig(file_save, dpi=300)

# %% >>> Save file to csv
file_save = path_probes / f'{name_probes}_ALL__with_blast_FLAPY__coverage.csv'
probes_summary.to_csv(file_save, sep=',', index=False)   
print(f'\n\nSUMMARY saved as {file_save}')

# %%
