#! /usr/bin/env python3

import os
import sys
from argparse import ArgumentParser

from Bio import SeqIO

parser = ArgumentParser(description='Generate a lot of FASTAs using list file')
parser.add_argument('-l', type=str, help='List file')
parser.add_argument('-f', type=str, nargs='+', help='FASTA file(s)')
parser.add_argument('-o', type=str, help='Output basename. Default `cluster`',
                    default='cluster')
parser.add_argument('-d', type=str, default='groups',
                    help='Output directory to be created')
parser.add_argument('--min', type=int, default=30,
                    help='Minimal cluster size')
parser.add_argument('--max', type=int, default=200,
                    help='Maximal cluster size')
parser.add_argument('--min_species', type=int, default=10,
                    help='Minimum species count in cluster')
args = parser.parse_args()

clusters = []
diatom_filehandles = []
counter = 0
os.mkdir(args.d)
# Assumes list to consist of tab-separated ID lists, one cluster in each line
for line in open(args.l):
    ids = line.rstrip().split('\t')
    species = len({x.split('|')[0] for x in ids})
    if args.min <= len(ids) <= args.max and species >= args.min_species:
        print(species)
        clusters.append(set(ids))
        diatom_filehandles.append(open('{0}/{1}{2}.diatoms.fasta'.format(args.d,
                                                                     args.o,
                                                                     counter),
                                       mode='w+'))
        counter += 1
print('Loaded {} clusters'.format(counter), flush=True, file=sys.stderr)

for fasta_file in args.f:
    fasta_iterator = SeqIO.parse(fasta_file, 'fasta')
    pwd = os.getcwd()
    os.chdir(args.d)
    for record in fasta_iterator:
        for index, cluster in enumerate(clusters):
            if record.id in cluster:
                SeqIO.write(record, diatom_filehandles[index], 'fasta')
    os.chdir(pwd)
print('Loaded diatom FASTAs', flush=True, file=sys.stderr)
