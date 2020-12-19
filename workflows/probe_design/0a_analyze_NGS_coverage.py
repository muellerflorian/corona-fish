
# THIS SCRIPT ALLOWS TO ANALYZE THE NGS GENOME COVERAGE.
# This has to be done only once, and it will produce a results file that can be directly used.
# 1. load all coverage files
# 2. Renormalize each coverage file by median
# 3. Calculate median of normalized coverage
# 4. Save this as `coverage.csv`. This file can then be used to select probes on highest coverage.  

# Script also produces a number of plots for inspection.

# %% Imports
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import json

# %% Read coverage files
path_coverage = Path(Path.cwd() / '..' / '..' / 'data' / 'coverage').resolve()
path_coverage_file = path_coverage / 'files'
# Read
cov_all_list = []
cov_info = {}
for idx, file_coverage in enumerate(path_coverage.glob('*.cov')):
    # Get name
    name_cov = f'cov_file__{idx}'
    print(f'Reading coverage file {file_coverage}, referenced as {name_cov}')
    # Read file
    cov_loop = pd.read_csv(file_coverage, sep='\t', header=None, names=['genome','pos', name_cov])
    cov_loop.drop(columns=['genome'], inplace=True)
    cov_all_list.append(cov_loop)
    # Same info about file
    cov_info[name_cov] = file_coverage.name

# Save file info to json    
file_save = path_coverage / 'cov_file_info.json'
with open(file_save, "w") as fh:
    json.dump(cov_info, fh, indent=4, sort_keys=True)  
    
# Generate one dataframe
cov_joined = cov_all_list[0]
cov_all_list = cov_all_list[1:]
for cov_add in cov_all_list:
    cov_joined = pd.merge(left=cov_joined, right=cov_add, how='outer', left_on='pos', right_on='pos')

# Renormalize each column (except position) by median
cov_joined_norm = cov_joined.div(cov_joined.median(axis=0), axis='columns') 
cov_joined_norm['pos'] = cov_joined['pos']

# >>>Get median coverage
cov_median_series = cov_joined_norm.drop(columns=['pos']).median(axis=1)
frame = {'pos':  cov_joined_norm['pos'], 'cov': cov_median_series }
cov_median = pd.DataFrame(frame) 
# Calculate rolling median
cov_median['cov_roll'] = cov_median['cov'].rolling(window=30).median()

file_save = path_coverage / 'coverage.csv' 
cov_median.to_csv(file_save, sep=',', index=False)   
print(f'\n\nSUMMARY saved as {file_save}')

# %% >>>>>>>>>>>>>>>>>>>> PLOTS FOR INSPECTION <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# %% Plot coverage over position
#  This is a bit slow since >25K data points are plotted

# Median normalized coverage
#sns.lineplot(x='pos', y='cov',data=cov_median)
sns.lineplot(x='pos', y='value', hue='variable', data=pd.melt(cov_median, ['pos']))
file_save = path_coverage / 'cov_norm_median.png'
plt.savefig(file_save, dpi=300)
plt.close()

# Raw data
sns.lineplot(x='pos', y='value', hue='variable', data=pd.melt(cov_joined, ['pos']))
file_save = path_coverage / 'cov_raw_pos.png'
plt.savefig(file_save, dpi=300)
plt.close()

# Renormalized data
sns.lineplot(x='pos', y='value', hue='variable', data=pd.melt(cov_joined_norm, ['pos']))
file_save = path_coverage / 'cov_norm_median_pos.png'
plt.savefig(file_save, dpi=300)
plt.close()

# %% Correlation between each columns

# PLot correlation matrix
cov_plot = cov_joined.drop(columns=['pos'])
cov_norm_plot = cov_joined_norm.drop(columns=['pos'])

f = plt.figure(figsize=(19, 15))
plt.matshow(cov_plot.corr(method ='pearson') , fignum=f.number)
plt.xticks(range(cov_plot.shape[1]), cov_plot.columns, fontsize=14, rotation=45)
plt.yticks(range(cov_plot.shape[1]), cov_plot.columns, fontsize=14)
cb = plt.colorbar()
cb.ax.tick_params(labelsize=14)
file_save = path_coverage / 'cov_raw_corr.png'
plt.savefig(file_save, dpi=300)
plt.close()

f = plt.figure(figsize=(19, 15))
plt.matshow(cov_norm_plot.corr(method ='pearson') , fignum=f.number)
plt.xticks(range(cov_norm_plot.shape[1]), cov_norm_plot.columns, fontsize=14, rotation=45)
plt.yticks(range(cov_norm_plot.shape[1]), cov_norm_plot.columns, fontsize=14)
cb = plt.colorbar()
cb.ax.tick_params(labelsize=14)
file_save = path_coverage / 'cov_norm_corr.png'
plt.savefig(file_save, dpi=300)
plt.close()