'''
Introduction
============

Module to collect various function design smFISH probes

INSTALLATION
^^^^^^^^^^^^

Blast ... requires command line tools

see https://biopython.org/DIST/docs/tutorial/Tutorial.html  Search for "Running BLAST locally"
and https://www.ncbi.nlm.nih.gov/books/NBK279690/


Usage
^^^^^
# Set up instance of probeDesigner class
probes = probeDesigner.smFISHprobes()
file_fasta ='./Probes_HIVcodingSequence--pol/Probes_HIVcodingSequence--pol_FILT_summary.fasta'

#### Perform blast - against RNA
file_xml     = './Probes_HIVcodingSequence--pol/results_blast__RNA.xml'
probes.performBlastnShort(file_xml,file_fasta,'rna')

#### Perform blast - against contig
file_xml     = './Probes_HIVcodingSequence--pol/results_blast__CONTIG.xml'
probes.performBlastnShort(file_xml,file_fasta,'all_contig')

### Summarize all of blast
probes.summarizeBlast()

### Save as csv file
fname   = './Probes_HIVcodingSequence--pol/results_blast__summary2.csv'
probes.writeSummaryCsv(fname)

'''

# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------

from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
import csv
from operator import itemgetter

# ---------------------------------------------------------------------------
# Exports
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Globals
# ---------------------------------------------------------------------------

__version__ = '0.0.1'
__author__ = 'Florian MUELLER'
__email___ = 'muellerf.research@gmail.com'


class smFISHprobes():
    '''Base class for AnnotationImporter object.'''

    def __init__(self):

        self.blast_summary_all = []
        self.blast_summary = []
        self.file_names_fasta = []
        self.file_names_xml = []
        self.blast_command = []
        self.blast = []

    def performBlastnShort(self, file_name_xml, file_name_fasta, db, evalue=5,max_target_seqs = 100,max_hsps = 1):
        '''
        Perform a blast with parameters for short sequences
        db ... locally stored sequence database (see blast+ for details)
            'allcontig_and_rna'  ... for blast against mouse genomic + transcript
            'rna'                ... for blast against human transcript
            'all_contig'         ... for blast against human genomic
        '''
        # Different outputs can be created; outfm = 7 gives text file, outfm = 5
        # gives XML
        # All options for command line: https://www.ncbi.nlm.nih.gov/books/NBK279675/
        blastn_cline = NcbiblastnCommandline(query=file_name_fasta,
                                             db=db,
                                             evalue=evalue,
                                             outfmt=5,
                                             out=file_name_xml,
                                             task='blastn-short',
                                             max_target_seqs = max_target_seqs,
                                             max_hsps = max_hsps)

        # Perform blast
        stdout, stderr = blastn_cline()

        # Summarize blast
        blast_summary = self.summarize_XML(file_name_xml)

        # Store variables
        self.blast_summary.append([blast_summary])
        self.blast_command.append([blastn_cline])
        self.file_names_fasta.append([file_name_fasta])
        self.file_names_xml.append([file_name_xml])
        

    def summarize_XML(self, file_xml, flag_print=False):
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

    def summarizeBlast(self):
        '''
        Summarize the different blast searches and order them by Query ID, and
        then length of alignment
        '''
        blast_summary_all = []

        for blast_res in self.blast_summary:

            # Fuse the two files and sort
            blast_summary_all.extend(blast_res[0])

        # Sort files
        blast_summary_all = sorted(
            blast_summary_all, key=itemgetter(3), reverse=True)
        blast_summary_all = sorted(blast_summary_all, key=itemgetter(0))

        self.blast_summary_all = blast_summary_all

    def writeSummaryCsv(self, file_name):
        ''' WRITE summarized results of all blast to a csv file '''

        # Open file with 'with' to guarantee that it will be closed
        with open(file_name, 'w') as fhandle:

            # Specify csv writer, csv.QUOTE_ALL - Quote everything, regardless
            # of type.
            writer = csv.writer(fhandle, quoting=csv.QUOTE_ALL)

            # Write information (header and alignment report)
            writer.writerow(('query_id', 'query_title', 'query_length', 'hsp_identities',
                             'hsp_gaps', 'category', 'alignment_title', 'hsp_score', 'hsp_expect'))
            last_query = ''
            for b in self.blast_summary_all:

                if len(last_query) == 0:
                    last_query = b[0]
                else:
                    if not(b[0] == last_query):
                        writer.writerow('')
                        last_query = b[0]
                writer.writerow(b)