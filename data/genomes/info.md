# Databases are  are created with blast suite. 
See `probe_design\_probe-design-overview.md` for more details.

Following commands have to be executed in the respective subfolders

## Beta-corona
```
makeblastdb -in beta-corona.fasta -dbtype nucl -title beta-corona
```

## Virus genomes
```
#!/bin/sh
makeblastdb -in chlamydophila-pneumoniae.fasta -dbtype nucl -title chlamydophila-pneumoniae -max_file_sz 500000 -parse_seqids
makeblastdb -in hpiv1.fasta -dbtype nucl -title hpiv1 -max_file_sz 500000 -parse_seqids
makeblastdb -in hpiv2.fasta -dbtype nucl -title hpiv2 -max_file_sz 500000 -parse_seqids
makeblastdb -in hpiv3.fasta -dbtype nucl -title hpiv3 -max_file_sz 500000 -parse_seqids
makeblastdb -in hpiv4.fasta -dbtype nucl -title hpiv4 -max_file_sz 500000 -parse_seqids
makeblastdb -in human-metapneumovirus.fasta -dbtype nucl -title human-metapneumovirus -max_file_sz 500000 -parse_seqids
makeblastdb -in mycoplasma-pneumoniae.fasta -dbtype nucl -title mycoplasma-pneumoniae -max_file_sz 500000 -parse_seqids
makeblastdb -in respiratory-syncytial-virus.fasta -dbtype nucl -title respiratory-syncytial-virus -max_file_sz 500000 -parse_seqids
makeblastdb -in tuberculosis.fasta -dbtype nucl -title tuberculosis -max_file_sz 500000 -parse_seqids
```

TODO: build for H3N2.fasta
