# Next steps

## Other output format
https://www.reddit.com/r/bioinformatics/comments/4ef5p8/how_to_filter_blast_results_using_biopython/

will only work with a local installation

So if you output tabular BLAST results (-outfmt 6), you can use Pandas to quickly
read, sort, filter your BLAST results without having to worry too much about
for loops, iterators and parsing properly. Pandas DataFrame is similar to the
R data.frame - basically a data table with a lot of nice convenience functions
for sorting, filtering, reading, writing, etc.

For example with tabular output file blast-out.b6:


## Make blast over the internet (no local installation required)
Starting point: https://www.biostars.org/p/129932/

``` python
from Bio.Blast import NCBIWWW
result_handle = NCBIWWW.qblast("blastn", "nt", some_sequence)
```
