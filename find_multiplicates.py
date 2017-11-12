#! /usr/bin/env python3

from argparse import ArgumentParser
from sys import stderr

from Bio import SeqIO

from phylome.blast_parser import parse_blast_file_to_hits, iterate_by_query
from phylome.multiplicates import is_duplicate


def duplicate_percentage(hit_iterable, overlap_cutoff, len_cutoff):
    """
    Return a percentage of duplicate hits in a list
    :param hit_iterable:
    :return:
    """
    duplicate_count = 0
    for hit in hit_iterable:
        if is_duplicate(hit, overlap_cutoff=overlap_cutoff,
                        len_cutoff=len_cutoff):
            duplicate_count += 1
    return duplicate_count/len(hit_iterable)


parser = ArgumentParser(description='Detect multiplicates based on BLAST hits')
parser.add_argument('-b', type=str, help='BLAST TSV file')
parser.add_argument('-f', type=str, default='',
                    help='FASTA file of queries')
parser.add_argument('-o', type=float, default=0.5,
                    help='Overlap cutoff. Default 0.5')
parser.add_argument('-l', type=int, default=50,
                    help='Length cutoff. Default 50')
parser.add_argument('-c', type=float, default=0.5,
                    help='Duplicate hit percentage. Default 0.5')
parser.add_argument('-n', type=str, default='multiples',
                    help='Output filename')
parser.add_argument('-s', action='store_true',
        help='Print species statistics. Assumes FASTA headers to be species|ID')
args = parser.parse_args()

if not args.b:
    print('Need a BLAST file', file=stderr)
    quit()
    
duplicates = set()
for query_group in iterate_by_query(parse_blast_file_to_hits(filename=args.b)):
    if duplicate_percentage(query_group, overlap_cutoff=args.o,
                            len_cutoff=args.l) >= args.c:
        duplicates.add(query_group[0].query_id)
print('Found {} duplicates'.format(len(duplicates)), file=stderr, flush=True)
if duplicates:
    species_mult = {}
    with open('{}.duplicate_list'.format(args.n), mode='w') as duplicate_list:
        for a in duplicates:
            print(a, file=duplicate_list)
            if args.s:
                species = a.split('|')[0]
                try:
                    species_mult[species] += 1
                except KeyError:
                    species_mult[species] = 1
    print('Generating FASTA...', file=stderr, flush=True)
    species_total = {}
    if args.f:
        with open('{}.duplicates.fasta'.format(args.n), mode='w') as outfasta:
            for record in SeqIO.parse(args.f, 'fasta'):
                if record.id in duplicates:
                    SeqIO.write(record, outfasta, 'fasta')
                if args.s:
                    species = record.id.split('|')[0]
                    try:
                        species_total[species] += 1
                    except KeyError:
                        species_total[species] = 1

if args.s:
    with open('{}.statistics.tsv'.format(args.n), mode='w') as stat_file:
        for species in sorted(species_total.keys()):
            if species not in species_mult.keys():
                species_mult[species] = 0
            print('\t'.join([species,
                     str(species_total[species]),
                     str(species_mult[species]),
                     str(species_mult[species]/species_total[species]*100)]),
                  file=stat_file)
