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
makeblastdb -in H1N1.fasta -dbtype nucl -title H1N1 -max_file_sz 500000 -parse_seqids
makeblastdb -in H3N2.fasta -dbtype nucl -title H3N2 -max_file_sz 500000 -parse_seqids
makeblastdb -in influenza_B.fasta -dbtype nucl -title influenza_B -max_file_sz 500000 -parse_seqids
makeblastdb -in influenza_C.fasta -dbtype nucl -title influenza_C -max_file_sz 500000 -parse_seqids
makeblastdb -in influenza_D.fasta -dbtype nucl -title influenza_D -max_file_sz 500000 -parse_seqids
makeblastdb -in rhenovirus_B.fasta -dbtype nucl -title rhenovirus_B -max_file_sz 500000 -parse_seqids
