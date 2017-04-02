#! /usr/bin/env python3

from argparse import ArgumentParser

parser = ArgumentParser(description='Cluster diatom sequences based on BLAST')
parser.add_argument('-b', type=str, help='BLAST tabular file')
parser.add_argument('-f', type=str, help='FASTA file')
args = parser.parse_args()

pass
