#! /usr/bin/env/python3

from argparse import ArgumentParser
from glob import glob
import os

parser = ArgumentParser('Create IQtree call lists')
parser.add_argument('-d', type=str, help='Alignment directory')
args = parser.parse_args()

command_mask = '/home/amorozov/tools/iqtree-1.6.1-Linux/bin/iqtree -s {0} -m LG+F+G -st AA -nt 1 >{0}.iqtree.stdout 2>{0}.iqtree.stderr'
if not os.path.isdir(args.d):
    raise ValueError('Incorrect path')
alignments = glob(args.d+'/*.aln.fasta')
for alignment in alignments:
    print(command_mask.format(os.path.abspath(alignment)))

