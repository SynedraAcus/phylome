#! /usr/bin/env python3

from argparse import ArgumentParser

import pymysql

from phylome.blast_parser import parse_blast_file_to_hits, iterate_by_query
from phylome.taxonomy import get_supertaxon_from_list

parser = ArgumentParser('Return best Archaeplastida hit')
parser.add_argument('-f', type=str, help='BLAST TSV file')
parser.add_argument('-d', type=str, help='mySQL DB name')
parser.add_argument('-o', type=str, help='mySQL host')
parser.add_argument('-u', type=str, help='mySQL username')
parser.add_argument('-p', type=str, help='mySQL password')
args = parser.parse_args()

cnx = pymysql.connect(user=args.u, host=args.o, password=args.p,
                                  database=args.d)
cursor = cnx.cursor()
taxid_request = 'select * from acc2taxid where accession in ({});'

for hit_list in iterate_by_query(parse_blast_file_to_hits(filename=args.f)):
    best = {'red': [], 'green': [], 'none': []}
    # Loading the relevant piece of acc2taxid
    r = taxid_request.format(', '.join(['"'+x.hit_id.split('.')[0]+'"' for x in hit_list]))
    cursor.execute(r)
    result = cursor.fetchall()
    acc2taxid = {}
    for x in result:
        acc2taxid[x[0]] = x[1]
    for hit in sorted(hit_list,
                      key=lambda x: min([hsp.evalue for hsp in x.hsps])):
        supertaxon = get_supertaxon_from_list(acc2taxid[hit.hit_id.split('.')[0]],
                                              [2166, 33090],
                                              cursor)
        if supertaxon:
            print(hit.hit_id, supertaxon)
            if supertaxon == 2166:
                best['red'].append(hit.hit_id)
            elif supertaxon == 33090:
                best['green'].append(hit.hit_id)
            break
        else:
            best['none'].append(hit.hit_id)
with open(args.f+'.reds', mode='w') as fh:
    for x in best['red']:
        print(x, file=fh)
with open(args.f+'.greens', mode='w') as fh:
    for x in best['green']:
        print(x, file=fh)
with open(args.f+'.none', mode='w') as fh:
    for x in best['none']:
        print(x, file=fh)
