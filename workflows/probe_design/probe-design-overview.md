# Probe-design for COV-2

- [Probe-design for COV-2](#probe-design-for-cov-2)
  - [Overview](#overview)
  - [Probe-design workflow](#probe-design-workflow)
    - [1. Create fasta file [Python]](#1-create-fasta-file-python)
      - [Workflow](#workflow)
    - [2. Design probes [R]](#2-design-probes-r)
      - [Workflow](#workflow-1)
    - [3. Blast sequences: LOCAL databases](#3-blast-sequences-local-databases)
      - [Build local database from downloaded fasta sequences](#build-local-database-from-downloaded-fasta-sequences)
      - [Workflow](#workflow-2)
    - [4.Blast sequence against human transcriptome](#4blast-sequence-against-human-transcriptome)
    - [5. Combine blast results](#5-combine-blast-results)
      - [Workflow](#workflow-3)
    - [6. Query probes](#6-query-probes)
      - [Workflow](#workflow-4)
  - [Additional information](#additional-information)
    - [Blast](#blast)
      - [Windows](#windows)
      - [BLASTn output format 6](#blastn-output-format-6)
      - [Local blast: Out of memory error: windows](#local-blast-out-of-memory-error-windows)
      - [Other output format for local blast](#other-output-format-for-local-blast)
      - [Make blast over the internet (no local installation required)](#make-blast-over-the-internet-no-local-installation-required)

## Overview

__Target__: **SARS-CoV-2** (NC_045512.2)
    Sequence of viruse stored in `data\Wuhan-Hu-1_2019.gb`

Target regions are define in the file `data\Wuhan-Hu-1_2019__target_regions.tsv`. This can be the entire genome, or sub-regions.

__Probe should be specific over other Coronaviruses:__

| beta-coronavirus | ID          | Included |
|------------------|-------------|----------|
| SARS             | NC_004718.3 | [x]      |
| MERS             | NC_019843.3 | [x]      |
| HKU1             | NC_006577.2 | [x]      |
| OC43             | NC_006213.1 | [x]      |
| NL63             | JX504050.1  | [x]      |
| 229E             | NC_002645.1 | [x]      |

`data\genomes\beta-corona`: genomes are stored in multi-fasta file, and blast databse build for it. 

__Additional specifity requirements:__

| Virus/pathogen                                             | ID          | Included |
|------------------------------------------------------------|-------------|----------|
| Mycobacterium tuberculosis                                 | NC_000962.3 | [x]      |
| Human parainfluenza virus type 1                           | NC_003461   | [x]      |
| Human parainfluenza virus type 2                           | NC_003443.1 | [x]      |
| Human parainfluenza virus type 3                           | NC_001796   | [x]      |
| Human parainfluenza virus type 4                           | NC_021928.1 | [x]      |
| Influenza A H3N2          | taxid 335341          | []       | 
| Influenza A H1N1 (such as A/Alaska/58/2017(H1N1))          |             | []       |
| Influenza B (such as B/Yamagata/16/88 and B/Victoria/2/87) |             | []       |
| Influenza C                                                |             | []       |
| Influenza D                                                |             | []       |
| Respiratory syncytial virus                                | NC_001803            | []       |
| Human metapneumovirus                                      | NC_039199             | []       |
| Rhinovirus/enterovirus                                     |             | []       |
| Mycoplasma pneumoniae                                      | NZ_CP010546            | []       |
| Chlamydophila pneumoniae                                   | NC_005043.1            | []       |

__Probes are aligned against different host organisms__:

| host                 | teaxid | tested |
|----------------------|--------|--------|
| Home sapiens         | 9606   | [x]    |
| African Green monkey | 60711  | []     |
| Mus musculus         | 10090  | []     |
| Hamsters             | 10026  | []     |
| Ferret               | 9669   | []     |


## Probe-design workflow

### 1. Create fasta file [Python]

Takes the file with the target regions, extracts their sequences, and
creates a separate fasta file for each target region type. If several sequences are 
defined for a target regions type, the script will create a multi-fasta file, e.g. 
multiple fasta entries in one file,

__Note__: fasta file names, sequence ids or description should **NOT** contain underscores. In Oligostan, underscores are used to split string to extract some meta-data.

#### Workflow

__Function__: `probe_design\1_create_fasta_of_target_region.py`

__Input__:

1. `data\Wuhan-Hu-1_2019.gb`
2. `data\Wuhan-Hu-1_2019_potential_target_regions.tsv`
   Different region types can be defined as well if needed (specified by the name in the first column). Currently supported are
     - `cov-2`
     - `spike-glycoprotein`
     - `target-region`
     - `primer-region`

  For a new region type, several target ranges (specified by their start and end position) can be specified.

  If you define a new region type, the Python script `1_create_fasta_of_target_region.py` to take new entry types into considerations.

__Output__:

1. `data\Wuhan-Hu-1_2019.gb`
2. `data\Wuhan-Hu-1_2019_potential_target_regions.tsv`

### 2. Design probes [R]

We use Oligostan for [probe-design](https://www.ncbi.nlm.nih.gov/pubmed/27599845).

See also the provided document `workflows\probe_design\Oligostan_documentation.doc` for more details. 

__IMPORTANT__: in the probe design script you have to update the path to the fasta sequence (line 8).

- Runs in **R** as a script, which designs probes against the fasta sequences created in step 1.
- Requires installation of some extra packages: `install.packages(c("ade4", "seqinr", "zoo"))`
- Takes **fasta files** as an input.
  - Multi-fasta is supported.
  - File-names, sequence ids or description can NOT contain underscores
  - When specifying path names under Windows, `\` has to be replaced by `/`

#### Workflow

__Function__: `probe_design\oligostan.r`

__Input__:

1. `data\Wuhan-Hu-1_2019.gb`
2. `data\Wuhan-Hu-1_2019_potential_target_regions.tsv`

__Output__:
Will create a folder `data\fasta\Probes__annotation-name`, where `region-name` is the name of the
regions given in step 1, e.g. `cov-2` when the annotated region was the default region covering the entire genome.

Contains a fasta file and a summary text file for all probes identified (`ALL`), and the probes passing the specified filters (`FILT`).

### 3. Blast sequences: LOCAL databases

Perform local blast against several reference genomes.

__IMPORTANT:__ if you add a new target region type, you have to update the path pointing to the create probe regions in step 2.

- Local version of blast + databases being installed.
- Some more information
  - https://biopython.org/DIST/docs/tutorial/Tutorial.html  search for "Running blast locally"
  - Blast+ user manual: https://www.ncbi.nlm.nih.gov/books/NBK279690/
  - https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download

- Downloads:
  - Blast+ suite: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
  - Instructions for different operating systems are provided.
    Important; blast commands have to be in system path. https://www.ncbi.nlm.nih.gov/books/NBK279671/
  - Database for human transcripts: https://ftp.ncbi.nlm.nih.gov/blast/db/

- Local databases can be build for alignment against any fasta sequence.

#### Build local database from downloaded fasta sequences

FASTA files for different viruses/pathogens/beta-coronaviruses are provided.

__Important__: you have to buidl the databases on your local installation. 

The general command is `makeblastdb -in beta-corona.fasta -dbtype nucl -title beta-corona -max_file_sz 500000 -parse_seqids`

- The `-parse_seqids` option is required to keep the original sequence identifiers.
- The `-max_file_sz` option helps to avoid an error under windows?

Has to be performed for each of the genomes were a blast is desired.

Build can provoke error messages under Windows. See misc below.

#### Workflow

__Function__: `probe_design\3_blast_local.py`

__Input__:

1. Fasta files with probe sequences
2. Local blastn databases

__Output__:

1. Blast hit files. http://www.metagenomics.wiki/tools/blast/blastn-output-format-6

### 4.Blast sequence against human transcriptome

1. Blast website: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch
2. Drag multi-fasta file with all sequences
3. Define blast search of interest
4. Download results as hit file (csv). We performed blast against
   1. Human: `blastn_human_transcriptome.csv`

### 5. Combine blast results

Results of all blast searches are combined, and a summary files with details of probe design and the blast results created.

#### Workflow

__Function__: `probe_design\5_blast_combine.py`

__Input__:

1. Fasta files from local search.
2. Fasta files from web search.

__Output__:

1. csv file with all information of probes.
2. Some plots showing mismatch distribution.

### 6. Query probes

Results as pandas dataframe. Query probes

#### Workflow

__Function__: `probe_design\6_probes_query.py`

__Input__:

1. Summary csv files

__Output__:

1. Fasta file with queried sequence.
2. Plot with distribution of probes along sequence.

## Additional information

### Blast

#### Windows

Installed under `C:\Program Files\NCBI\blast-2.10.0+`

`Path` was automatically update with bin directory: `C:\Program Files\NCBI\blast-2.10.0+\bin`

#### BLASTn output format 6

Obtained when downloading blastn internet searches as a **hit table**.
BLASTn tabular output format 6.

- `qseqid`: query (e.g., unknown gene) sequence id
- `sseqid`: subject (e.g., reference genome) sequence id
- `pident`: percentage of identical matches
- `length`: alignment length (sequence overlap)
- `mismatch`: number of mismatches
- `gapopen`: number of gap openings
- `qstart`: start of alignment in query
- `qend`: end of alignment in query
- `sstart`: start of alignment in subject
- `send`: end of alignment in subject
- `evalue`: expect value
- `bitscore`: bit score

__Note__: additional fields could be added. From blastn -help:

`qseq` means Aligned part of query sequence
`sseq` means Aligned part of subject sequence

Note that qseq/sseq may contain gaps (`-` characters).

#### Local blast: Out of memory error: windows

This can create an error on windows (https://github.com/flu-crew/octoFLU/issues/16). Can be resolved by 
by setting an environmental variable `BLASTDB_LMDB_MAP_SIZE=1000000`

Steps to create or modify environment variables are summarized below:

1. Search for "Environment Variables" and Select "Edit environmental ...".
2. Click the "New" button under the "User variable for ..." panel
3. Type the environment variable `BLASTDB_LMDB_MAP_SIZE` and set its value to `1000000`
4. Click "OK" to close the prompts

#### Other output format for local blast

https://www.reddit.com/r/bioinformatics/comments/4ef5p8/how_to_filter_blast_results_using_biopython/

will only work with a local installation

So if you output tabular BLAST results (-outfmt 6), you can use Pandas to quickly
read, sort, filter your BLAST results without having to worry too much about
for loops, iterators and parsing properly. Pandas DataFrame is similar to the
R data.frame - basically a data table with a lot of nice convenience functions
for sorting, filtering, reading, writing, etc.

For example with tabular output file blast-out.b6:

#### Make blast over the internet (no local installation required)

Starting point: https://www.biostars.org/p/129932/

``` python
from Bio.Blast import NCBIWWW
result_handle = NCBIWWW.qblast("blastn", "nt", some_sequence)
```
