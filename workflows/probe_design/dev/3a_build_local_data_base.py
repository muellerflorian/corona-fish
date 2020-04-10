
# %%
import os
import Path
import subprocess

list_files = subprocess.run([f"cd {path_viruses}"])
# %%  Path to genomes

path_genomes = Path(Path.cwd() / '..' / '..' / 'data' / 'genomes').resolve()
path_viruses = str(path_genomes / 'viruses')

out = os.system(f"cd {path_viruses}")
print(out)
# %%

for fasta_genome in path_viruses.glob(f'*.fasta'):
    print(f'Buidling genome: {fasta_genome}')

        


# %%


    # >> Web blast
    blast_all = []
    for file_blast in path_blast.glob(f'*.{ext}'):
        print(f'Reading blast results: {file_blast}')
        file_name = file_blast.name