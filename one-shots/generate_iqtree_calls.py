#! /usr/bin/env python3

import os
from argparse import ArgumentParser
from glob import glob

parser = ArgumentParser('Create IQtree call lists')
parser.add_argument('-d', type=str, help='Alignment directory')
args = parser.parse_args()

command_mask = '/home/amorozov/tools/iqtree-1.6.1-Linux/bin/iqtree -s {0} -m LG+F+G -st AA -nt 1 >{0}.iqtree.stdout 2>{0}.iqtree.stderr'
if not os.path.isdir(args.d):
    raise ValueError('Incorrect path')
alignments = glob(args.d+'/*.aln.fasta')
trees = [x.split('.')[0]+'.dataset.fasta.aln.fasta' for x in glob(args.d+'/*.treefile')]
to_call = filter(lambda x: x not in trees, alignments)
for alignment in to_call:
    print(command_mask.format(os.path.abspath(alignment)))

