
# TODO: include missing blasts (Mycobacterium tuberculosis, Influenza A)
# TODO: important or not to include genetic variability within virus?

# %% Imports
from pathlib import Path
import importlib
from covfish import probe_designer
import sys

# %% Blast against local databases
path_probes = Path(Path.cwd() / '..' / '..' / 'data' / 'fasta' / 'Probes__cov-2').resolve()
path_genomes = Path(Path.cwd() / '..' / '..' / 'data' / 'genomes').resolve()

importlib.reload(probe_designer)
file_probes = path_probes / 'Probes__cov-2_ALL.fasta'

# Just to make sure that the target sequences are correct
db = path_genomes / 'cov-2' / 'cov-2.fasta'
probe_designer.blast_probes(file_probes, db, min_identity=0.75,  add_new_line=False)

#  Blast against other beta-coronaviruses
db = path_genomes / 'beta-corona' / 'beta-corona.fasta'
probe_designer.blast_probes(file_probes, db, min_identity=0.75,  add_new_line=False)

# >>> Blast against other viruses
path_db_blast = path_genomes / 'viruses'
blast_summary_list = []
for db in path_db_blast.glob('*.fasta'):
    print(f'Performing local blast against: {db}')
    try:
        blast_loop = probe_designer.blast_probes(file_probes, db, min_identity=0.75, add_new_line=False)
        blast_summary_list.append([blast_loop])
    except:
        print('Blast failed. Maybe database is not build?')
        e = sys.exc_info()[0]
        print(e)


# %%
