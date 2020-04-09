
# %% Imports
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
import csv
from operator import itemgetter
from pathlib import Path

#%% Function definition to summarize results
def summarize_XML(file_xml, flag_print=False):
    ''' Summarize results of a blast search stored in a xml file '''
    blast_summary = []

    for record in NCBIXML.parse(open(file_xml)):

        # Add leading zeros to query ID
        query_number = record.query_id.split('_')[1].zfill(3)
        query_id_new = 'Query_'+query_number

        if flag_print is True:
            print()
            print(record.query_id)
            print(record.query)

        for alignment in record.alignments:
            for hsp in alignment.hsps:

                if flag_print is True:
                    print('')
                    print('****Alignment****')
                    print('sequence:', alignment.title)
                    print('length:', alignment.length)
                    print('score:', hsp.score)
                    print('identities:', hsp.identities)
                    print('gaps:', hsp.gaps)
                    print('e-value:', hsp.expect)
                    print(hsp.query)
                    print(hsp.match)
                    print(hsp.sbjct)

                # Analyse alignment title
                #    https://en.wikipedia.org/wiki/RefSeq
                #    http://www.nslc.wustl.edu/elgin/genomics/bio4342/1archives/2006/AccReference.pdf
                category = '-'

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


                blast_loop = (query_id_new, record.query, record.query_length,
                                hsp.identities, hsp.gaps, category, alignment.title, hsp.score, hsp.expect)
                blast_summary.append(blast_loop)

    return blast_summary
def write_summary_csv(file_name, blast_summary, min_identity = 0.7, add_new_line=False):
    """[summary]
    
    Parameters
    ----------
    file_name : [type]
        [description]
    blast_summary : [type]
        [description]
    min_identity : float, optional
        [description], by default 0.7
    """

    # Open file with 'with' to guarantee that it will be closed
    with open(file_name, 'w', newline='') as fhandle:

        # Specify csv writer, csv.QUOTE_ALL - Quote everything, regardless of type.
        writer = csv.writer(fhandle, quoting=csv.QUOTE_ALL)

        # Write information (header and alignment report)
        writer.writerow(('query_id', 'query_title', 'query_length', 'hsp_identities',
                            'hsp_gaps', 'category', 'alignment_title', 'hsp_score', 'hsp_expect'))
        current_query = ''
        
        # Read row by row, add empty row after new target sequence
        for b in blast_summary:
            new_line_to_add = False
            if add_new_line:
                if not (b[0] == current_query):
                    writer.writerow('')
                    new_line_to_add = True
                    current_query = b[0]
                
            if b[3]/b[2] > min_identity:
                if new_line_to_add:
                    writer.writerow('')
                writer.writerow(b)
#%%

# Blast list of probes against local db
def blast_probes(file_fasta, db, min_identity=0.7, add_new_line=False):

    path_save = file_fasta.parents[0] / 'blast'
    if not path_save.is_dir():
        path_save.mkdir()
        
    # Define function call
    db_name = db.stem
    file_xml = path_save / f'{db_name}.xml'
    file_csv = path_save / f'{db_name}.csv' 
    
    blastn_cline = NcbiblastnCommandline(query=file_fasta,
                                        db= db,
                                        evalue=5,
                                        outfmt=5,
                                        out=file_xml,
                                        task='blastn-short',
                                        max_target_seqs = 100,
                                        max_hsps = 1)
    # Perform blast ans summarize XML file
    stdout, stderr = blastn_cline()

    blast_summary = summarize_XML(file_xml=file_xml, flag_print=False)
    write_summary_csv(file_csv, blast_summary, min_identity = min_identity, add_new_line=add_new_line )

#%% Blast against local databases

file_fasta=Path(r'D:\Work\Documents\Projects\collaborations\corona-fish\probe-design\oligostan\length_26-32\Probes__cov-2\Probes__cov-2_ALL_summary.fasta')

db_blast = [Path(r'D:\Work\GitHub\projects\corona-fish\data\genomes\cov-2\cov-2.fasta'),
            Path(r'D:\Work\GitHub\projects\corona-fish\data\genomes\corona_family\corona_family.fasta'),
            Path(r'D:\Work\GitHub\projects\corona-fish\data\genomes\tuberculosis\tuberculosis.fasta')]


# Loop over databases
for db in db_blast:
    blast_probes(file_fasta, db, min_identity=0.7,  add_new_line=False)
    




# %% READ xml file
file_xml = r'D:\Work\Documents\Projects\collaborations\corona-fish\probe-design\oligostan\length_22-26\Probes_target-reg-seq\Probes_target-reg-seq_ALL_blast.xml'
blast_summary = summarize_XML(file_xml=file_xml,
               flag_print=False)


# %% Write to csv
file_save = r'D:\Work\Documents\Projects\collaborations\corona-fish\probe-design\oligostan\length_22-26\Probes_target-reg-seq\Probes_target-reg-seq_ALL_blast_summary.csv'
write_summary_csv(file_save, blast_summary, min_identity = 0.65)

# %%
