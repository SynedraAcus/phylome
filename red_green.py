#! /usr/bin/env python3

import mysql.connector
from argparse import ArgumentParser
from phylome.blast_parser import parse_blast_file_to_hits
from phylome.taxonomy import get_supertaxon_from_list

parser = ArgumentParser(description='Read BLAST file, count red&green hits')
parser.add_argument('-b', type=str, help='BLAST file')
args = parser.parse_args()

cnx = mysql.connector.connect(user='amorozov', host='10.12.0.251', password='a279#8y9z', database='amorozov')
cursor = cnx.cursor()

taxa = {}
for hit in parse_blast_file_to_hits(args.b):
    if not hit.query_id in taxa.keys():
        taxa[hit.query_id] = [0, 0]
    hit_id = hit.hit_id.split('.')[0]
    cursor.execute('SELECT ncbi_taxon_id FROM acc2taxid WHERE accession=\'{}\';'
                   .format(hit_id))
    mysql_answer = cursor.fetchall()
    if not mysql_answer:
        # Silently passing 
        continue
    taxon_id = mysql_answer[0][0]
    supertaxon = get_supertaxon_from_list(taxon_id, [2763, 15432], cursor)
    if supertaxon:
        if supertaxon == 2763:
            taxa[hit.query_id][0] += 1
        else:
            taxa[hit.query_id][1] += 1
for taxon in taxa.keys():
    print('\t'.join((taxon, str(taxa[taxon][0]), str(taxa[taxon][1]))))
