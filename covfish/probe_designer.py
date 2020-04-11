
# Imports
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
import csv
from operator import itemgetter
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import json

# Perform local blast
def blast_probes(file_fasta, db, path_save, add_new_line=False, max_target_seqs=100, result_fields=None):
    """[summary]

    Parameters
    ----------
    file_fasta : [type]
        [description]
    db : [type]
        [description]
    folder_sub : [type]
        [description]
    min_identity : float, optional
        [description], by default 0.7
    add_new_line : bool, optional
        [description], by default False
    max_target_seqs : int, optional
    result_fields : str with additional format specifiers
        Additional format specifiers that should be added to hit list file.
        http://www.metagenomics.wiki/tools/blast/blastn-output-format-6

    Returns
    -------
    [type]
        [description]
    """
    print(f'Performing local blast against: {db}')

    # Path to save results
    if not path_save.is_dir():
        path_save.mkdir(parents=True)

    # File to save results
    db_name = db.stem    
    file_save = path_save / f'{db_name}.txt'

    # Output format
    if result_fields:
        outfmt = f'6 {result_fields}'
    else:
        outfmt = '6'

    # Perform blast
    blastn_cline = NcbiblastnCommandline(query=file_fasta,
                                        db=db,
                                        strand='both',
                                        evalue=5,
                                        outfmt=outfmt,
                                        out=file_save,
                                        task='blastn-short',
                                        max_target_seqs=max_target_seqs,
                                        max_hsps=1)
    # Perform blast ans summarize XML file
    stdout, stderr = blastn_cline()

# Read blast results and summarize them
def blast_summarize(path_blast, file_identifier, sep, ext, path_results, refseq_keep=None):
    """[summary]

    Parameters
    ----------
    path_blast : pathlib object
        Path to folder containg blast search results.
    file_identifier : pathlib object
        Full file name of json file with identifiers for blast search
    sep : str
        Column separator in blast hit file
    ext: str
        File extension of blast hit file
    path_results : pathlib object
        Where to store results
    refseq_keep : list of str
        Refseq entries that should be kept

    Returns
    -------
    [list]
        List of panda dataframes with best blast hits.
    """

    # Column names of blast hit file
    names_col = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'hsend', 'evalue', 'bitscore']

    # Check if folder exists
    if not path_blast.is_dir():
        print(f'Folder with blast results not found {path_blast}')

    # Open file-idenifiers for blast searches
    with open(file_identifier, 'r') as fp:
        blast_name_identifier = json.load(fp)

    # >> Web blast
    blast_all = []
    for file_blast in path_blast.glob(f'*.{ext}'):
        print(f'Reading blast results: {file_blast}')
        file_name = file_blast.name

        # Internet blast against human  (separator is ')
        blast_results = pd.read_csv(file_blast, sep=sep, header=None, names=names_col)

        # Keep only transcripts
        if refseq_keep:
            blast_results = blast_results[blast_results.sseqid.str.contains('|'.join(refseq_keep))]

        # Keep only top hits
        df_sorted = blast_results.sort_values(by='bitscore')
        best_blast_hits = df_sorted.drop_duplicates('qseqid', keep='last').sort_index(0)

        # Rename columns
        if not (file_name in blast_name_identifier):
            ident = file_blast.stem
            print(f'No identifier for result file  {file_name} found. Will use file stem {ident} instead.')

        else:
            ident = blast_name_identifier[file_name]

        # Plot results of mismatches
        if len(best_blast_hits) > 0:
            file_save = path_results / f'blast_{ident}_mismatch.png'
            hexplot = sns.jointplot("length", "mismatch", data=best_blast_hits, kind="hex")       
            plt.subplots_adjust(left=0.2, right=0.8, top=0.8, bottom=0.2)  # shrink fig so cbar is visible
            cbar_ax = hexplot.fig.add_axes([.85, .25, .05, .4])  # x, y, width, height
            plt.colorbar(cax=cbar_ax)
            plt.savefig(file_save, dpi=300)
            plt.close()
        else:
            print('No blast hits for this search')

        # Rename columns
        best_blast_hits.rename(columns={'sseqid': f'{ident}_sseqid',
                                        'pident': f'{ident}_pident',
                                        'length': f'{ident}_length',
                                        'mismatch': f'{ident}_mismatch',
                                        'gapopen': f'{ident}_gapopen',
                                        'qstart': f'{ident}_qstart',
                                        'qend': f'{ident}_qend',
                                        'sstart': f'{ident}_sstart',
                                        'hsend': f'{ident}_hsend',
                                        'evalue': f'{ident}_evalue',
                                        'bitscore': f'{ident}_bitscore'
                                        },
                               inplace=True)

        # Add results
        blast_all.append(best_blast_hits)

    return blast_all


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