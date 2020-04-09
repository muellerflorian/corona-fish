
# TODO: include missing blasts (Mycobacterium tuberculosis, Influenza A)
# TODO: important or not to include genetic variability within virus?

# %% Imports
from pathlib import Path
import importlib
from covfish import probe_designer

# %% Blast against local databases
path_probes = Path(Path.cwd() / '..' / '..' / 'data' / 'fasta' / 'Probes__cov-2').resolve()
path_genomes = Path(Path.cwd() / '..' / '..' / 'data' / 'genomes').resolve()

importlib.reload(probe_designer)
file_fasta = path_probes / 'Probes__cov-2_ALL_summary.fasta'

# Just to make sure that the target sequences are correct
db = path_genomes / 'cov-2' / 'cov-2.fasta'
probe_designer.blast_probes(file_fasta, db, min_identity=0.75,  add_new_line=False)

# Blast against other databases
db_blast = [path_genomes / 'beta-corona' / 'beta-corona.fasta',
            path_genomes / 'tuberculosis' / 'tuberculosis.fasta']

blast_summary_list = []
for db in db_blast:
    blast_loop = probe_designer.blast_probes(file_fasta, db, min_identity=0.75, add_new_line=False)
    blast_summary_list.append([blast_loop])


# %%
