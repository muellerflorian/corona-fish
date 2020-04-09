
# %% Sequence Input
#  It is for use when the handle contains one and only one record, which is returned as a single SeqRecord object. 
#  If there are no records, or more than one, then an exception is raised:

from Bio import SeqIO
record = SeqIO.read("single.fasta", "fasta")


# %% Sequence manipulations

# Reverse complement
my_dna.reverse_complement()


# %% Sequence Output
# Use the function Bio.SeqIO.write(...), which takes a complete set of SeqRecord objects (either as a list, or an iterator), an output file handle (or in recent versions of Biopython an output filename as a string) and of course the file format:
from Bio import SeqIO
sequences = ...  # add code here
with open("example.fasta", "w") as output_handle:
    SeqIO.write(sequences, output_handle, "fasta")
    
    
