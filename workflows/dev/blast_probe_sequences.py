# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %%
## REQUIREMENT
#  - BioPython
#  - local version of blast + databases being installed: 
#        - https://biopython.org/DIST/docs/tutorial/Tutorial.html  search for "Running blast locally"
#        - https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download
#        - Need to be added to system path to be found in command line


# IMPORTS 
import os
import sys

# Change this path to the folder containing the python code
sys.path.append('/Volumes/PILON_HD2/fmueller/Documents/code/development/fishquant/scr')
import probeDesigner


# %%
# Define fasta file
path_fasta = '/Volumes/PILON_HD2/fmueller/Documents/Projects/smFISH_collaborations/Allain/Probes/SHANK3/EX21-22'
file_fasta ='Probes_Shank3-21-22_FILT_summary.fasta'

# Set up instance of probeDesigner class
probes = probeDesigner.smFISHprobes()


# %%
# Perform blast -  data-base have to be downloaded
#  All options for command line: https://www.ncbi.nlm.nih.gov/books/NBK279675
#     'allcontig_and_rna'  ... for blast against mouse genomic + transcript  
#     'rna'                ... for blast against human transcript
#     'all_contig'         ... for blast against human genomic


# %%
# HUMAN BLAST - RNA and genomic contigs. Results will be be fused
probes.performBlastnShort(os.path.join(path_fasta,'results_blast__RNA.xml'),os.path.join(path_fasta,file_fasta),'rna',max_target_seqs = 20,max_hsps=1)
probes.performBlastnShort(os.path.join(path_fasta,'results_blast__CONTIG.xml'),os.path.join(path_fasta,file_fasta),'all_contig',max_target_seqs = 20,max_hsps=1)


# %%
# MOUSEBLAST
probes.performBlastnShort(os.path.join(path_fasta,'results_blast.xml'),os.path.join(path_fasta,file_fasta),'allcontig_and_rna',max_target_seqs = 20,max_hsps=1)


# %%
# Summarize and save
probes.summarizeBlast()
probes.writeSummaryCsv(os.path.join(path_fasta,'results_blast__summary.csv'))



