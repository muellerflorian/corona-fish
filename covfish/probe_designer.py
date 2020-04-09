
# Imports
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
import csv
from operator import itemgetter
from pathlib import Path

# Function definition to summarize results
def summarize_XML(file_xml, flag_print=False):
    ''' 
    Summarize results of a blast search stored in a xml file 
    
    High-scoring Segment Pair (HSP) class: https://biopython.org/DIST/docs/api/Bio.Blast.Record.HSP-class.html
    '''
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
                category = 'NA'

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


                blast_loop = (query_id_new, record.query, record.query_length, hsp.align_length, 
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
        writer.writerow(('query_id', 'query_title', 'query_length', 'alignment_length', 'hsp_identities',
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


# Blast list of probes against local db
def blast_probes(file_fasta, db, min_identity=0.7, add_new_line=False):

    path_save = file_fasta.parents[0] / 'blast'
    if not path_save.is_dir():
        path_save.mkdir()
        
    # Define function call
    db_name = db.stem
    
    # For xml format
    #outfmt = 5 
    #file_save = path_save / f'{db_name}.xml'
    
    # Tabular format
    outfmt = 6
    file_save =  path_save / f'{db_name}.txt'
    
    blastn_cline = NcbiblastnCommandline(query=file_fasta,
                                        db= db,
                                        evalue=5,
                                        outfmt=outfmt,
                                        out=file_save,
                                        task='blastn-short',
                                        max_target_seqs = 100,
                                        max_hsps = 1)
    # Perform blast ans summarize XML file
    stdout, stderr = blastn_cline()

    #file_csv = path_save / f'{db_name}.csv' 
    #blast_summary = summarize_XML(file_xml=file_xml, flag_print=False)
    #write_summary_csv(file_csv, blast_summary, min_identity=min_identity, add_new_line=add_new_line)
    #return blast_summary

def blast_results_combine(blast_summary):
        '''
        Summarize the different blast searches and order them by Query ID, and
        then length of alignment
        '''
        blast_summary_all = []

        for blast_res in blast_summary:

            # Fuse the two files and sort
            blast_summary_all.extend(blast_res[0])

        # Sort files
        blast_summary_all = sorted(blast_summary_all, key=itemgetter(3), reverse=True)
        blast_summary_all = sorted(blast_summary_all, key=itemgetter(0))

        return blast_summary_all