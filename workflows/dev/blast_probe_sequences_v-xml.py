
# %% Imports
import sys
from pathlib import Path
import pandas as pd 
import importlib

sys.path.insert(0, str((Path.cwd() / '..' / 'src').resolve()) )
import probe_designer

#%% Blast against local databases
importlib.reload(probe_designer)
file_fasta = Path(r'D:\Work\Documents\Projects\collaborations\corona-fish\probe-design\oligostan\length_26-32\Probes__cov-2\Probes__cov-2_ALL_summary.fasta')

# Just to make sure that the target sequences are correct
db = Path(r'D:\Work\GitHub\projects\corona-fish\data\genomes\cov-2\cov-2.fasta')
probe_designer.blast_probes(file_fasta, db, min_identity=0.75,  add_new_line=False)

# Blast against other databases 
db_blast = [ Path(r'D:\Work\GitHub\projects\corona-fish\data\genomes\corona_family\corona_family.fasta'),
             Path(r'D:\Work\GitHub\projects\corona-fish\data\genomes\tuberculosis\tuberculosis.fasta')]

blast_summary_list = []
for db in db_blast:
    blast_loop = probe_designer.blast_probes(file_fasta, db, min_identity=0.75, add_new_line=False)
    blast_summary_list.append([blast_loop])
    
blast_summary = probe_designer.blast_results_combine(blast_summary_list)
path_save = file_fasta.parents[0] / 'blast'
if not path_save.is_dir():
   path_save.mkdir()

file_csv = path_save / '_blast_off_target.csv'
probe_designer.write_summary_csv(file_csv, blast_summary, min_identity=0.75, add_new_line=False)

# %% Add  blast results to probe list

# >>>> DF with blast result
df_blast = pd.DataFrame(blast_summary, columns = ['query_id', 'query_title', 'query_length', 'alignment_length', 'hsp_identities',
                            'hsp_gaps', 'category', 'alignment_title', 'hsp_score', 'hsp_expect'])

df_blast['perc_align_corona_fam'] = 100 * (df_blast['hsp_identities'] / df_blast['query_length'])
#df_blast['n_mismatch'] = df_blast['query_length'] - df_blast['hsp_identities']

aggregation_functions = {'perc_align_corona_fam': 'max'}
df_blast_agg = df_blast.groupby(df_blast['query_title']).aggregate(aggregation_functions)

# >>>> DF with list of probes
file_summary = Path(r'D:\Work\Documents\Projects\collaborations\corona-fish\probe-design\oligostan\length_26-32\Probes__cov-2\Probes__cov-2_ALL.txt')
probes_summary = pd.read_csv(file_summary, sep='\t')

# Create matching column name with blast results
probes_summary['query_title'] = probes_summary['ProbesNames'] + '--' + probes_summary['theStartPos'].astype('str') + '--' + probes_summary['theEndPos'].astype('str')

# >>> Join data DFs
df_joined = pd.merge(left=probes_summary, right=df_blast_agg, how='left', left_on='query_title', right_on='query_title')

# Save as csv
file_save = Path(r'D:\Work\Documents\Projects\collaborations\corona-fish\probe-design\oligostan\length_26-32\Probes__cov-2\Probes__cov-2_ALL_with_blast.csv')
df_joined.to_csv(file_save, sep=',')


# %% Blast against human genome and transcriptome
# csv file. individual queries are not separated.
file_blast_human = Path(r'D:\Work\Documents\Projects\collaborations\corona-fish\probe-design\oligostan\length_26-32\Probes__cov-2\blast\blastn_human_transcriptome.csv')
blastn_human = pd.read_csv(file_blast_human, sep=',', header=0,
                           names = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore'])



# %% Could also be read Biopython. This uses a generator object
# Here hit list has to exported as text file. This text file contains comments indicated with a #
from Bio import SearchIO
file_blast_human = Path(r'D:\Work\Documents\Projects\collaborations\corona-fish\probe-design\oligostan\length_26-32\Probes__cov-2\blast\blastn_human_transcriptome.txt')
qresult = SearchIO.parse(file_blast_human, 'blast-tab', comments=True)


#%%

if 'NM_' in alignment.title:
    category = 'mRNA'

elif 'XM_' in alignment.title:
    category = 'mRNA_pred'

elif 'NR_' in alignment.title:
    category = 'ncRNA'

elif 'XR' in alignment.title:
    category = 'ncRNA_pred'

elif 'NT' in alignment.title:
    category = 'genomic'

elif 'NW' in alignment.title:
    category = 'genomic-automated'

# %% Impport results from internet blast 

# READ xml file
file_xml = r'D:\Work\Documents\Projects\collaborations\corona-fish\probe-design\oligostan\length_22-26\Probes_target-reg-seq\Probes_target-reg-seq_ALL_blast.xml'
blast_summary = summarize_XML(file_xml=file_xml,
               flag_print=False)

file_save = r'D:\Work\Documents\Projects\collaborations\corona-fish\probe-design\oligostan\length_22-26\Probes_target-reg-seq\Probes_target-reg-seq_ALL_blast_summary.csv'
write_summary_csv(file_save, blast_summary, min_identity = 0.65)

# %%
