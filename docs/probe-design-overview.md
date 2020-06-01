# Probe-design for SARS-Cov-2

- [Probe-design for SARS-Cov-2](#probe-design-for-sars-cov-2)
  - [Overview](#overview)
    - [Target **SARS-CoV-2** (NC_045512.2)](#target-sars-cov-2-nc_0455122)
    - [NGS coverage](#ngs-coverage)
      - [Background](#background)
    - [Other beta-coronaviruses](#other-beta-coronaviruses)
    - [Other viruses/pathogens](#other-virusespathogens)
    - [Probes are (web) blasted against possible host genomes](#probes-are-web-blasted-against-possible-host-genomes)
  - [Probe-design workflow](#probe-design-workflow)
    - [1. Create fasta file for target sequence(s) [Python]](#1-create-fasta-file-for-target-sequences-python)
      - [Workflow](#workflow)
    - [2. Design probes [R]](#2-design-probes-r)
      - [Workflow](#workflow-1)
    - [3. Blast sequences: LOCAL databases](#3-blast-sequences-local-databases)
      - [Local blast install](#local-blast-install)
      - [Build local database from downloaded fasta sequences](#build-local-database-from-downloaded-fasta-sequences)
      - [Local blast: Out of memory error: windows](#local-blast-out-of-memory-error-windows)
      - [Workflow](#workflow-2)
    - [4.Blast sequence against human transcriptome](#4blast-sequence-against-human-transcriptome)
    - [5. Combine blast results](#5-combine-blast-results)
      - [Workflow](#workflow-3)
    - [6. NGS coverage of probes](#6-ngs-coverage-of-probes)
      - [Workflow](#workflow-4)
    - [7. Query probes](#7-query-probes)
      - [Workflow](#workflow-5)
  - [Additional information](#additional-information)
    - [Blast](#blast)
      - [Windows](#windows)
      - [BLASTn output format 6](#blastn-output-format-6)

## Overview

### Target **SARS-CoV-2** (NC_045512.2)

- **Genome** sequence is stored here: `data\Wuhan-Hu-1_2019.gb`
- For a **blast** search against this genome, we use: `data\genomes\cov2\cov2.fasta`
- For a **blast** against >2500 **aligned cov-2 genomes**, we use: `data\genomes\cov2\cov2_aligned.fasta`. Please see `data\genomes\cov2\gisaid_acknowledgements.txt` for the acknowledgments.
- **Target regions** are define in the file `data\Wuhan-Hu-1_2019__target_regions.tsv`. This can be the entire genome, or sub-regions.

### NGS coverage

Probes that target regions for having highest NGS coverage are selected. All relevant data is in the folder `data\coverage`.

If new coverage files are added, they first have to be analyzed with the script `workflows\probe_design\0a_analyze_NGS_coverage.py`.

#### Background

> Coverage refers to the average number of times a single base is read during a sequencing run.

What is Coverage in NGS? In theory, a sequencing run will geenrate reads that sample a genome randomly and independently. 
However, these reads are not distributed equally across an entire genome; some bases are covered by fewer, some by more reads 
than the average coverage.

As a simplification, coverage can be assumed to be **proportional to the amount of RNA species in the sample**,
and should hence allow to pinpoint towards regions more likely to be present in a sample (as genomes/subgenomes/transcripts).

### Other beta-coronaviruses

Identified probes should not be specific against any other beta-coronavirus.

`data\genomes\beta-corona`: genomes are stored in multi-fasta file.

| beta-coronavirus | ID          | Included |
|------------------|-------------|----------|
| SARS             | NC_004718.3 | [x]      |
| MERS             | NC_019843.3 | [x]      |
| HKU1             | NC_006577.2 | [x]      |
| OC43             | NC_006213.1 | [x]      |
| NL63             | JX504050.1  | [x]      |
| 229E             | NC_002645.1 | [x]      |

### Other viruses/pathogens

Identified probes should not be specific against other pathogens/viruses with similar symptomes.

`data\genomes\viruses`: genomes are stored in separate fasta file.

| Virus/pathogen                                             | ID               | Included |
|------------------------------------------------------------|------------------|----------|
| Mycobacterium tuberculosis                                 | NC_000962.3      | [x]      |
| Human parainfluenza virus type 1                           | NC_003461        | [x]      |
| Human parainfluenza virus type 2                           | NC_003443.1      | [x]      |
| Human parainfluenza virus type 3                           | NC_001796        | [x]      |
| Human parainfluenza virus type 4                           | NC_021928.1      | [x]      |
| Respiratory syncytial virus                                | NC_001803        | [x]      |
| Human metapneumovirus                                      | NC_039199        | [x]      |
| Mycoplasma pneumoniae                                      | NZ_CP010546      | [x]      |
| Chlamydophila pneumoniae                                   | NC_005043.1      | [x]      |
| Influenza A H3N2                                           |                  | [x]      |
| Influenza A H1N1 (such as A/Alaska/58/2017(H1N1))          | [txid2043069](https://www.ncbi.nlm.nih.gov/nuccore/?term=txid2043069[Organism:noexp])      | [x]      |
| Influenza B (such as B/Yamagata/16/88 and B/Victoria/2/87) | [txid416674](https://www.ncbi.nlm.nih.gov/nuccore/?term=txid416674[Organism:noexp])       | [x]      |
| Influenza C                                                | [GCF_000856665.10](https://www.ncbi.nlm.nih.gov/nuccore/1250175392,1250175391,1250175390,1250175389,1250175388,211910015,52630357) | [x]      |
| Influenza D                                                | [GCF_002867775.1](https://www.ncbi.nlm.nih.gov/nuccore/1328406502,1328406500,1328406498,1328406496,1328406494,1328406492,1328406490)  | [x]      |
| Rhinovirus/enterovirus                                     | NC_038312.1      | [x]      |
  
### Probes are (web) blasted against possible host genomes

Identified probes should not be specific against transcriptome of commonly used host organisms.

Only tanscripts are kept from blast hits, e.g. Refseq annotations starting with : 'NM_', 'XM_', 'NR', 'XR_'

| Host                 | Identifier | Tested |
|----------------------|------------|--------|
| Home sapiens         | Human G+T  | [x]    |
| Mus musculus         | Mouse G+T  | [x]    |
| African Green monkey | 60711      | [x]    |
| Hamsters             | 10026      | [x]    |
| Ferret               | 9669       | [x]    |

## Probe-design workflow

### 1. Create fasta file for target sequence(s) [Python]

Takes the file with the target regions, extracts their sequences, and
creates a separate fasta file for each target region type. If several sequences are 
defined for a target regions type, the script will create a multi-fasta file, e.g. 
multiple fasta entries in one file.

In this function, you can define if the fasta sequence should be the **reverse complement**.

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

  If you define a new region type, you have to modify the the Python script `1_create_fasta_of_target_region.py` to take 
  new entry types into considerations.

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

Perform local blast against several reference genomes. Results are stored as blast hit files with output format 6
(http://www.metagenomics.wiki/tools/blast/blastn-output-format-6).

Blast is performed with the option `blastn-short`, which is specifically optimized for shorter sequences (<50 bases). This changes the following
parameters

- `word-size`: 2
- `gapopen`: 5
- `gapextend`: 2
- `reward`: 1
- `penalty`: -3
  
__IMPORTANT__: you have to build the databases on your local installation (see below).

#### Local blast install

You can download the blast+ suite from this [link](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/). 
Instructions for different operating systems are provided.

__Important__: blast commands have to be in system path as described [here](https://www.ncbi.nlm.nih.gov/books/NBK279671/).

- Local databases can be build for alignment against any fasta sequence.

#### Build local database from downloaded fasta sequences

FASTA files for different viruses/pathogens/beta-coronaviruses are provided. blast data-base has to be build
for each new local installation.

Build commands are listed in `data\genomes\info.md`

The **general command** is `makeblastdb -in beta-corona.fasta -dbtype nucl -title beta-corona -max_file_sz 500000 -parse_seqids`

- The `-parse_seqids` option is required to keep the original sequence identifiers.
- The `-max_file_sz` option helps to avoid an error under windows?

Has to be performed for each of the genomes where a blast is desired.

Build can provoke error messages under Windows. See  below.

#### Local blast: Out of memory error: windows

This can create an error on windows (https://github.com/flu-crew/octoFLU/issues/16). Can be resolved by 
by setting an environmental variable `BLASTDB_LMDB_MAP_SIZE=1000000`

Steps to create or modify environment variables are summarized below:

1. Search for "Environment Variables" and Select "Edit environmental ...".
2. Click the "New" button under the "User variable for ..." panel
3. Type the environment variable `BLASTDB_LMDB_MAP_SIZE` and set its value to `1000000`
4. Click "OK" to close the prompts

#### Workflow

__Function__: `probe_design\3_blast_local.py`

__Input__:

1. Fasta files with probe sequences: `data\fasta\Probes__cov-2\Probes__cov-2_ALL.fasta`
2. Local blastn databases: `data\genomes\`

__Output__:

1. Blast hit files for each data-basesstored in `\data\fasta\Probes__cov-2\blast\local`

### 4.Blast sequence against human transcriptome

1. [Blast website](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch)
2. **Enter query sequence**: Upload file with probe sequences**: `data\fasta\Probes__cov-2\Probes__cov-2_ALL.fasta`
3. **Choose Search Set**:
   1. `Database`: `Standard databases (nr etc.)`
   2. `Organism`: add taxid of host organism
4. **Program Selection**: `Somewhat similar sequences`
5. Perform Blast.
6. **Download results** as hit file (csv). We performed blast against

### 5. Combine blast results

Results of all blast searches are combined, and a summary files with details of probe design and the blast results created.

The file `data\blast_identifiers.json` permits to specify an identifier for a blast search. If absent, the file-name will be used.

#### Workflow

__Function__: `probe_design\5_blast_combine.py`

__Input__:

1. Fasta files from local search.
2. Fasta files from web search.

__Output__:

1. csv file with all information of probes.
2. Some plots showing mismatch distribution.

### 6. NGS coverage of probes

We calculate for each probe the NGS coverage (summarized with the script (`workflows\probe_design\0a_analyze_NGS_coverage.py`).

#### Workflow

__Function__: `workflows\probe_design\6_coverage.py`

__Input__:

1. csv file with all information of probes and blast results.

__Output__:

1. csv file with all information of probes, blast results, and probe coverage
2. Plots showing coverage of each probe.

### 7. Query probes

Results are provided as a csv file combining all information about probe design, blast analysis, and probe coverage. 

Probes are queried for

- GC content
- Passing sequence filters
- Highest NGS coverage
- Maximum overlap with different cov-2 genomes
- Minimum overlap with other beta-corona viruses, other viruses/pathogens causing similar symptomes.

#### Workflow

__Function__: `probe_design\7_probes_query.py`

__Input__:

1. Summary csv file.

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

`qseq` stands for  aligned part of **q**uery sequence
`sseq` stands for aligned part of **s**ubject sequence

Note that qseq/sseq may contain gaps (`-` characters).
