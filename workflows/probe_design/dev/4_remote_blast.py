# ATTEMPT for remote blast, but this seems very slow.


# %% Some ressources
#
# >>> qblast : blast over the internet (no local installation required)
#
# INPUT
#   https://www.biostars.org/p/249818/
#   - qblast will accept fasta or a single sequence. 
#        > fastafile = open('sequence.fna')
#        > result_handle = NCBIWWW.qblast('blastn', 'nt', fastafile.read())
#        > fastafile.close()
#  
# 
# Blast searching:  http://www.dalkescientific.com/writings/NBN/blast_searching.html
#
#
#

# Returns: The qblast function returns a cStringIO.StringIO instance.
# - Python class which acts like a file 
# - Instead of getting the data from the file system it gets the data from a string. 
# - If you want to see the result you can read from this file (or "file-like") instance using it's read() method
# 
# Parser


#### Make b
# #Starting point: https://www.biostars.org/p/129932/

# %% Imports
from pathlib import Path
import importlib
from Bio.Blast import NCBIWWW

# %% Blast against local databases
path_probes = Path(Path.cwd() / '..' / '..' / 'data' / 'fasta' / 'Probes__cov-2').resolve()
file_probes = path_probes / 'Probes__cov-2_ALL_summary.fasta'

# %% Perform blastn against transcriptomes

# >>> Homo sapiens
fastafile = open(file_probes)
result_handle = NCBIWWW.qblast('blastn', 'refseq_rna', fastafile.read())
fastafile.close()


# %%
#``` python

parser = NCBIWWW.BlastParser()
result_handle = NCBIWWW.qblast("blastn", "refseq_rna", "GAGTCACAGTTTGCTGTTTCTTCTGTCTCT",format_type = 'Text')
result = parser.parse(qblast_output) # http://www.bioinformatics.org/bradstuff/bp/tut/Tutorial003.html

#refseq_rna as a database ... 

#https://www.biostars.org/p/58457/
#NCBIWWW.qblast(program="blastn", database="refseq_genomic", sequence=input_sequence, entrez_query="txid7227[ORGN]")

with open("XML_Files/48_50.xml", "w") as save_to:
    save_to.write(result_handle.read())
    result_handle.close()