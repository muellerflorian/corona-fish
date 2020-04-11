# Databases are  are created with blast suite. 
See `probe_design\_probe-design-overview.md` for more details.

Following commands have to be executed in the respective subfolders

## Beta-corona
```
makeblastdb -in beta_corona.fasta -dbtype nucl -title beta_corona -max_file_sz 500000 -parse_seqids
```

## Cov-2
__Note__: aligned cov-2 contains 2531 cov-2 sequences.

```
makeblastdb -in cov2.fasta -dbtype nucl -title cov2 -max_file_sz 500000 -parse_seqids
makeblastdb -in cov2_aligned.fasta -dbtype nucl -title cov2_aligned -max_file_sz 500000 -parse_seqids
```

## Virus genomes
```
makeblastdb -in chlamydophila_pneumoniae.fasta -dbtype nucl -title chlamydophila_pneumoniae -max_file_sz 500000 -parse_seqids
makeblastdb -in hpiv1.fasta -dbtype nucl -title hpiv1 -max_file_sz 500000 -parse_seqids
makeblastdb -in hpiv2.fasta -dbtype nucl -title hpiv2 -max_file_sz 500000 -parse_seqids
makeblastdb -in hpiv3.fasta -dbtype nucl -title hpiv3 -max_file_sz 500000 -parse_seqids
makeblastdb -in hpiv4.fasta -dbtype nucl -title hpiv4 -max_file_sz 500000 -parse_seqids
makeblastdb -in human_metapneumovirus.fasta -dbtype nucl -title human_metapneumovirus -max_file_sz 500000 -parse_seqids
makeblastdb -in mycoplasma_pneumoniae.fasta -dbtype nucl -title mycoplasma_pneumoniae -max_file_sz 500000 -parse_seqids
makeblastdb -in respiratory_syncytial_virus.fasta -dbtype nucl -title respiratory_syncytial_virus -max_file_sz 500000 -parse_seqids
makeblastdb -in tuberculosis.fasta -dbtype nucl -title tuberculosis -max_file_sz 500000 -parse_seqids
makeblastdb -in H1N1.fasta -dbtype nucl -title H1N1 -max_file_sz 500000 -parse_seqids
makeblastdb -in H3N2.fasta -dbtype nucl -title H3N2 -max_file_sz 500000 -parse_seqids
makeblastdb -in influenza_B.fasta -dbtype nucl -title influenza_B -max_file_sz 500000 -parse_seqids
makeblastdb -in influenza_C.fasta -dbtype nucl -title influenza_C -max_file_sz 500000 -parse_seqids
makeblastdb -in influenza_D.fasta -dbtype nucl -title influenza_D -max_file_sz 500000 -parse_seqids
makeblastdb -in rhenovirus_B.fasta -dbtype nucl -title rhenovirus_B -max_file_sz 500000 -parse_seqids
```