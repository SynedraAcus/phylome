#! /usr/bin/env python3

from Bio import SeqIO
from argparse import ArgumentParser

parser = ArgumentParser(description='Filter FASTA file by the percentage of X and lowercase letters')
parser.add_argument('-f', type=str, help='FASTA file')
parser.add_argument('-n', type=float, help='Max acceptable proportion of X and lowercase')
args = parser.parse_args()

def lowercase_count(line):
    """
    A count of lowercase letters and/or X'es in a string 
    :param line: 
    :return: 
    """
    count = 0
    for letter in line:
        if letter.islower() or letter == 'X':
    count += 1

with open('{}.filtered'.format(args.f), mode='w') as outfasta:
    for record in SeqIO.parse(args.f, 'fasta'):
        if lowercase_count(record.seq)/len(record) <= args.n:
            SeqIO.write(record, outfasta, 'fasta')
