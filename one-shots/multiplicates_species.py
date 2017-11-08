#! /usr/bin/env python3

from argparse import ArgumentParser
from Bio import SeqIO

parser = ArgumentParser('Count sequences from different species')
parser.add_argument('-f', type=str, help='Sequence FASTA file')
parser.add_argument('-t', type=str, help='Sequence ID to taxon mapping')
args = parser.parse_args()

species = {}
species_totals = {}
for line in open(args.t):
    a = line.rstrip().split('\t')
    species[a[0]] = a[1]
    if a[1] in species_totals.keys():
        species_totals[a[1]] += 1
    else:
        species_totals[a[1]] = 1
counts = {x: 0 for x in species_totals.keys()}
for record in SeqIO.parse(args.f, format='fasta'):
    counts[species[record.id]] += 1
for a in counts.keys():
    print('\t'.join((a, str(counts[a]), str(counts[a]/species_totals[a]))))