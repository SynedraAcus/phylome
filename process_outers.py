#! /usr/bin/env python3

# process outers FASTA(s): for each supplied FASTA file, iterate over records,
# removing any diatoms (ncbi taxid 2836) and checking whether there are red or
# green (NCBI taxid 2763 and 33090, respectively). Archaeplastida are written to
# {original_fasta}.archaeplastida, others are written to {original_fasta}.rest
# The script also writes a log file as the TSV. The format is:
# original_fasta    #arche    #rest
#
# Also some spam in sys.stderr


from argparse import ArgumentParser


def process_fasta(filename, host='127.0.0.1', username='root',
                  password='password', database='biosql',
                  excluded = [2836], nonreduced=[2763, 33090]):
    """
    Process a single FASTA file. Returns a log line
    Log consists of filename, counts of all nonreduced (in the order in which
    they are in the kwarg) and the count of non-excluded, non-nonreduced seqs
    :param filename: 
    :return: 
    """
    # These are gonna be multiprocessed and each process will import
    # independently anyway.
    from Bio import SeqIO
    from phylome.taxonomy import get_supertaxon_from_list
    import mysql.connector
    records = {x.id.split('.')[0] for x in SeqIO.parse(filename, 'fasta')}
    cnx = mysql.connector.connect(user=username, host=host, password=password,
                                  database=database)
    cursor = cnx.cursor()
    # In a separate variable solely for the readability
    taxid_request = 'select * from acc2taxid where accession in ({});'
    r = taxid_request.format(', '.join(['\"'+x+'\"' for x in records.keys()]))
    cursor.execute(r)
    result = cursor.fetchall()
    tax2seq = {}
    for x in result:
        if x[1] in tax2seq:
            tax2seq[x[1]].append(x[0])
        else:
            tax2seq[x[1]] = [x[0]]
    nonreduced_seqid = []
    rest_seqid = []
    nonreduced_counts = {x: 0 for x in nonreduced}
    for taxid in tax2seq:
        supertaxon = get_supertaxon_from_list(taxid, excluded+nonreduced, cursor)
        if supertaxon:
            # The taxid is either from excluded or from nonreduced
            # excluded are silently passed
            if supertaxon in nonreduced:
                nonreduced_seqid += tax2seq[supertaxon]
                nonreduced_counts[supertaxon] += 1
            else:
                rest_seqid += tax2seq[supertaxon]
    SeqIO.write([records[x] for x in nonreduced_seqid],
                '{}.archaeplastida'.format(filename), 'fasta')
    SeqIO.write([records[x] for x in rest_seqid],
                '{}.rest'.format(filename), 'fasta')
    return '\t'.join([filename] +
                     [str(nonreduced_counts[x]) for x in nonreduced] +
                      [str(len(rest_seqid))])
            
            
            
    
    
    

parser = ArgumentParser('Process outers FASTA(s)')
parser.add_argument('-f', nargs='*', type=str,
                    help='FASTA file(s)')
parser.add_argument('-l', type=str, help='log file')
parser.add_argument('-h', type=str, help='MySQL host')
parser.add_argument('-u', type=str, help='MySQL username')
parser.add_argument('-p', type=str, help='MySQL password')
parser.add_argument('-t', type=int, help='Process number', default=2)
args = parser.parse_args()
