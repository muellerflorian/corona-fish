# %% Imports
from pathlib import Path
import csv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# %% Read cov-2 sequence from gb file
file_corona_seq = Path(Path.cwd() / '..' / '..' /'data' / 'Wuhan-Hu-1_2019.gb').resolve()
file_target_regions = Path(Path.cwd() / '..' / '..' / 'data' / 'Wuhan-Hu-1_2019__target_regions.tsv').resolve()
reverse_complement = False

# %% Read sequence file

# Corona sequence
if not file_corona_seq.is_file():
    print('File with corona sequence not found.')
    print(f'{file_corona_seq}')

corona_seq_record = SeqIO.read(file_corona_seq, "gb")
corona_seq = corona_seq_record.seq

# %% Analyze target sequences and create fasta files for each sequence type
# TODO: add option to create reverse complement to target coding sequences.

if not file_target_regions.is_file():
    print('File with target regions not found.')
    print(f'{file_target_regions}')

primer_seq = []
target_region_seq = []
spike_glyco_seq = []
cov_2_seq = []

counter_primer = 1
counter_region = 1


with open(file_target_regions) as tsvfile:
    target_regions = csv.DictReader(tsvfile, dialect='excel-tab')

    for target_region in target_regions:

        # Extract sequence
        # Note 1: Indexing in Python starts at 0
        # Note 2: Index range in python doesn't include last index, so no change here
        index_start = int(target_region['start']) - 1
        index_end = int(target_region['End'])

        seq_loop = corona_seq[index_start: index_end]
        seq_rec = SeqRecord(seq_loop)

        if reverse_complement:
            seq_rec = seq_rec.reverse_complement()

        seq_rec.description = f'start-{index_start+1}--end-{index_end}'

        if target_region['Name'] == 'primer-region':
            seq_rec.id = f'primer-{counter_primer}'
            primer_seq.append(seq_rec)
            counter_primer += 1

        elif target_region['Name'] == 'target-region':
            seq_rec.id = f'target-region-{counter_region}'
            target_region_seq.append(seq_rec)
            counter_region += 1

        elif target_region['Name'] == 'spike-glycoprotein':
            seq_rec.id = f'spike-glycoprotein'
            spike_glyco_seq.append(seq_rec)

        elif target_region['Name'] == 'cov-2':
            seq_rec.id = f'cov-2'
            cov_2_seq.append(seq_rec)

# %% Write target sequences to fasta file
path_save = Path(Path.cwd() / '..' / '..' / 'data'/ 'probes' ).resolve()
if not path_save.is_dir():
    path_save.mkdir(parents=True)
    print(f'Fasta sequences will be saved in folder {path_save}')

if reverse_complement:
    suffix = '--RevComp.fasta'
else:
    suffix = '.fasta'

if len(cov_2_seq) > 0:
    file_save = path_save / f'cov-2{suffix}'
    with open(file_save, "w") as output_handle:
        SeqIO.write(cov_2_seq, output_handle, "fasta")  

if len(primer_seq) > 0:
    file_save = path_save / f'primer-seq{suffix}'
    with open(file_save, "w") as output_handle:
        SeqIO.write(primer_seq, output_handle, "fasta")

if len(target_region_seq) > 0:
    file_save = path_save / f'target-reg-seq{suffix}' 
    with open(file_save, "w") as output_handle:
        SeqIO.write(target_region_seq, output_handle, "fasta")

if len(spike_glyco_seq) > 0:
    file_save = path_save / f'spike-glyco-seq{suffix}'    
    with open(file_save, "w") as output_handle:
        SeqIO.write(spike_glyco_seq, output_handle, "fasta")
 
# %%
