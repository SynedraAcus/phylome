#! /usr/bin/env python3

from Bio import SeqIO
from argparse import ArgumentParser

parser = ArgumentParser(description='Generate a lot of FASTAs using list file')
parser.add_argument('-l', type=str, help='List file')
parser.add_argument('-f', type=str, nargs='+', help='FASTA file(s)')
parser.add_argument('-o', type=str, help='Output basename. Default `cluster`',
                    default='cluster')
args = parser.parse_args()

clusters = []
filehandles = []
counter = 0
for line in open(args.l):
    counter += 1
    exec('clusters.append({})'.format(line.rstrip()))
    if not type(clusters[-1]) is set:
        raise ValueError('Incorrect line {}'.format(line))
    filehandles.append(open('{0}{1}.fasta'.format(args.o,
                                                  counter),
                            mode='w'))

records = [[] for x in clusters]
for fasta_file in args.f:
    for record in SeqIO.parse(fasta_file, 'fasta'):
        for index, cluster in enumerate(clusters):
            if record.id in cluster:
                records[index].append(cluster)

for index, record_group in enumerate(records):
    SeqIO.write(record_group, filehandles[index])

