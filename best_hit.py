#! /usr/bin/env python3

import sys
from argparse import ArgumentParser

import pymysql

from phylome.blast_parser import parse_blast_file_to_hits, iterate_by_query
from phylome.taxonomy import get_supertaxon_from_list


def cached_supertaxon(taxon, supertaxa, cursor, cache={}):
    """
    Get supertaxon from a list.
    Takes a taxon, a list of (potential) supertaxa and a pymysql cursor.
    If a taxon is a descendant of one of these supertaxa, returns that
    supertaxon, otherwise returns None.
    This function is cached and stores all the requests made during this
    program's run. Assumes that it wouldn't be called with different supertaxa
    lists during a single run.
    :return:
    """
    if taxon not in cache:
        cache[taxon] = get_supertaxon_from_list(taxon, supertaxa,
                                                cursor)
    return cache[taxon]
    
    
parser = ArgumentParser('Return best Archaeplastida hit')
parser.add_argument('-f', type=str, help='BLAST TSV file')
parser.add_argument('-d', type=str, help='mySQL DB name')
parser.add_argument('-o', type=str, help='mySQL host')
parser.add_argument('-u', type=str, help='mySQL username')
parser.add_argument('-p', type=str, help='mySQL password')
parser.add_argument('--best', action='store_true',
                help='Store only if the Archaeplastida hit is the very best')
args = parser.parse_args()

cnx = pymysql.connect(user=args.u, host=args.o, password=args.p,
                                  database=args.d)
cursor = cnx.cursor()
taxid_request = 'select * from acc2taxid where accession in ({});'

# Opening handles
red_handle = open(args.f + '.reds', mode='w')
green_handle = open(args.f + '.greens', mode='w')
none_handle = open(args.f + '.none', mode='w')
best = {'red': [], 'green': [], 'none': []}
for hit_list in iterate_by_query(parse_blast_file_to_hits(filename=args.f)):
    # Loading the relevant piece of acc2taxid
    r = taxid_request.format(', '.join(['"'+x.hit_id.split('.')[0]+'"' for x in hit_list]))
    cursor.execute(r)
    result = cursor.fetchall()
    acc2taxid = {}
    for x in result:
        acc2taxid[x[0]] = x[1]
    l = sorted(hit_list, key=lambda x: min([hsp.evalue for hsp in x.hsps]))
    if args.best:
        hit = l[0]
        try:
            taxon = acc2taxid[hit.hit_id.split('.')[0]]
        except KeyError:
            print('Unknown sequence ID {}'.format(hit.hit_id))
            continue
        try:
            supertaxon = cached_supertaxon(taxon, [2166, 33090], cursor)
        except ValueError:
            print('Unknown ID {}'.format(taxon), file=sys.stderr)
            continue
        if supertaxon:
            print(hit.query_id, supertaxon)
            if supertaxon == 2166:
                print(hit.query_id, file=red_handle)
                red_handle.flush()
                best['red'].append(hit.query_id)
            elif supertaxon == 33090:
                best['green'].append(hit.query_id)
                print(hit.query_id, file=green_handle)
                green_handle.flush()
        else:
            best['none'].append(hit.query_id)
            print(hit.query_id, file=none_handle)
            none_handle.flush()
    else:
        for hit in l:
            try:
                taxon = acc2taxid[hit.hit_id.split('.')[0]]
            except KeyError:
                print('Unknown sequence ID {}'.format(hit.hit_id))
                continue
            try:
                supertaxon = cached_supertaxon(taxon, [2166, 33090], cursor)
            except ValueError:
                print('Unknown ID {}'.format(taxon), file=sys.stderr)
                continue
            if supertaxon:
                print(hit.query_id, supertaxon)
                if supertaxon == 2166:
                    print(hit.query_id, file=red_handle)
                    red_handle.flush()
                    best['red'].append(hit.query_id)
                elif supertaxon == 33090:
                    best['green'].append(hit.query_id)
                    print(hit.query_id, file=green_handle)
                    green_handle.flush()
                break
            else:
                best['none'].append(hit.query_id)
                print(hit.query_id, file=none_handle)
                none_handle.flush()

print('Red\t{}\nGreen\t{}\nNone\t{}'.format(len(best['red']),
                                            len(best['green']),
                                            len(best['none'])))
for handle in (red_handle, green_handle, none_handle):
    handle.close()
