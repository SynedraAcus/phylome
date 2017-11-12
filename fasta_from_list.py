#! /usr/bin/env python3

import os
from argparse import ArgumentParser

from Bio import SeqIO

parser = ArgumentParser(description='Generate a lot of FASTAs using list file')
parser.add_argument('-l', type=str, help='List file')
parser.add_argument('-f', type=str, nargs='+', help='FASTA file(s)')
parser.add_argument('-o', type=str, help='Output basename. Default `cluster`',
                    default='cluster')
parser.add_argument('-d', type=str, default='groups',
                    help='Output directory to be created')
parser.add_argument('--min', type=int, default=50,
                    help='Minimal cluster size')
parser.add_argument('--max', type=int, default=200,
                    help='Maximal cluster size')
args = parser.parse_args()

clusters = []
filehandles = []
counter = 0
os.mkdir(args.d)
# Assumes list to consist of tab-separated ID lists, one cluster in each line
for line in open(args.l):
    ids = line.rstrip().split('\t')
    if args.min <= len(ids) <= args.max:
        clusters.append(set(ids))
        filehandles.append(open('{0}/{1}{2}.fasta'.format(args.d, args.o,
                                                      counter),
                           mode='w+'))
        counter += 1

for fasta_file in args.f:
    fasta_iterator = SeqIO.parse(fasta_file, 'fasta')
    pwd = os.getcwd()
    os.chdir(args.d)
    for record in fasta_iterator:
        for index, cluster in enumerate(clusters):
            if record.id in cluster:
                SeqIO.write(record, filehandles[index], 'fasta')
    os.chdir(pwd)
