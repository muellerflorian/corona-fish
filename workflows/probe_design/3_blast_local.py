
# %% Specify which data-set should be analyzed
name_probes = 'Probes__cov-2'  # Genomic probes
name_probes = 'Probes__cov-2--RevComp'  # Replication intermediate

# %% Imports
from pathlib import Path
import importlib
from covfish import probe_designer
import sys

# %% Blast against local databases
importlib.reload(probe_designer)

# Folders with probe sequences
path_data = Path(Path.cwd() / '..' / '..' / 'data' ).resolve()
path_probes = path_data / 'probes' / f'{name_probes}'
file_probes = path_probes / f'{name_probes}_ALL.fasta'  # Sequence without FLAP
file_probes_FLAP = path_probes / f'{name_probes}_ALL_FLAPY.fasta'  # Sequence without FLAP

# Folder with blast data-base and to store results
path_genomes = Path(Path.cwd() / '..' / '..' / 'data' / 'genomes').resolve()
path_save = path_probes / 'blast'

#  >>> Just to make sure that the target sequences are correct
db = path_genomes / 'cov2' / 'cov2.fasta'
probe_designer.blast_probes(file_probes, db, path_save=path_save / 'local_cov2', add_new_line=False)

#  >>> Blast against known cov-2 sequences
result_fields = 'qseqid sseqid qlen pident length mismatch gapopen qstart qend sstart send evalue bitscore'
db = path_genomes / 'cov2' / 'cov2_aligned.fasta'
n_genomes = db.read_text().count('>') 
probe_designer.blast_probes(file_probes, db,  path_save=path_save / 'local_cov2', add_new_line=False, max_target_seqs=n_genomes, result_fields=result_fields)

# >>>  Blast against other beta-coronaviruses
db = path_genomes / 'beta_corona' / 'beta_corona.fasta'
probe_designer.blast_probes(file_probes_FLAP, db,  path_save=path_save / 'local', add_new_line=False)

# >>> Blast against other viruses
path_db_blast = path_genomes / 'viruses'
blast_summary_list = []
for db in path_db_blast.glob('*.fasta'):
    try:
        blast_loop = probe_designer.blast_probes(file_probes_FLAP, db,  path_save=path_save / 'local', add_new_line=False)
        blast_summary_list.append([blast_loop])
    except:
        print('Blast failed. Maybe database is not build?')
        e = sys.exc_info()[0]
        print(e)


# %%
