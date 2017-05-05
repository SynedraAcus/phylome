#! /usr/bin/env python3

import re
import sys
from argparse import ArgumentParser


def yield_query_groups(line_iterable):
    """
    Yield GFF3 lines grouped by query.
    Assumes file to be query_sorted
    :param line_iterable: 
    :return: 
    """
    query = None
    r = []
    for line in line_iterable:
        a = line.split('\t')
        if not query:
            query = a[0]
        if query == a[0]:
            r.append(line)
        else:
            yield r
            query = a[0]
            r = [line]
    return r


def extract_group_gos(line_group, regex):
    """
    Extract all GOs.
    Accepts regex as a param to save time on regex compilation
    :param line_group: 
    :param regex: 
    :return: 
    """
    gos = set()
    for line in line_group:
        go_list = re.findall(regex, line)
        gos.update(go_list)
    return gos


def read_obo_file(obo_file):
    """
    Read an OBO file and return the dict of {GO: GO_name}
    :param obo_file: 
    :return: 
    """
    r = {}
    running_id = ''
    for line in open(obo_file):
        if 'id: ' in line:
            running_id = line.rstrip()[7:]
        if 'name: ' in line:
            r[running_id] = line.rstrip()[6:]
    return r

parser = ArgumentParser(description='Count GO terms in the GFF3 file and read'+\
                                    'their names from OBO')
parser.add_argument('-g', type=str, help='GFF3 file, stripped of comments and FASTA')
parser.add_argument('-o', type=str, help='OBO file')
args = parser.parse_args()
regex = re.compile('"GO\:(\d+)"')

print('Reading GFF3...', file=sys.stderr)
gos = {}
for line_group in yield_query_groups(open(args.g)):
    query_gos = extract_group_gos(line_group, regex)
    for go in query_gos:
        if go in gos:
            gos[go] += 1
        else:
            gos[go] = 1
print('Reading OBO...', file=sys.stderr)
go_names = read_obo_file(args.o)
for x in sorted(gos.keys(), key=lambda x: gos[x], reverse=True):
    try:
        print('{}\t{}'.format(go_names[x], gos[x]))
    except Keyerror:
        print('GO:{}\t{}'.format(x, gos[x]))
