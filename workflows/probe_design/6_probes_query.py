
# %% Imports
import pandas as pd
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# %% Query probes

path_probes = Path(Path.cwd() / '..' / '..' / 'data' / 'fasta' / 'Probes__cov-2').resolve()
file_results = path_probes / 'Probes__cov-2_ALL_with_blast.csv'
probes_summary_load = pd.read_csv(file_results)

# Replace NaN by zeros for easier analysis
probes_summary_load.h_length[probes_summary_load.h_length == 'NAN'] = 0.0
probes_summary_load.cf_length[probes_summary_load.cf_length == 'NAN'] = 0.0

probes_summary_load = probes_summary_load.astype({'cf_length': 'float64'})
probes_summary_load = probes_summary_load.astype({'h_length': 'float64'})

print(probes_summary_load.query(' h_length<20 ').shape[0])
print(probes_summary_load.query('GCFilter==1 & NbOfPNAS>2 & h_length<20 & cf_length<20').shape[0])


# %% Plot probe positions
df_query = probes_summary_load.query('GCFilter==1 & NbOfPNAS>2 & h_length<20 & cf_length<20')
file_save = path_probes / 'probe_positions.png'

plt.figure(figsize=(10,1))
plt.plot(df_query['theStartPos'].values, np.ones((df_query.shape[0],1)), '|', color='black')
plt.tight_layout()
plt.savefig(file_save, dpi =300)

# %%
