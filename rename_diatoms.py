#! /usr/bin/env python3

import shutil
from argparse import ArgumentParser
from tempfile import NamedTemporaryFile

from Bio import SeqIO


def rename_sequence(record, renamer=lambda x: x)
    record.description = renamer(record.description)
    record.id = ''
    record.name = ''
    #  return for readability only
    return record


def rename_all(file, renamer, in_place=False):
    with NamedTemporaryFile(mode='w+t') as outfile:
        for record in SeqIO.parse(file, format='fasta'):
            record = rename_sequence(record)
            SeqIO.write(record, outfile, 'fasta')
        if in_place:
            outname = file
        else:
            outname = '{}.renamed'.format(file)
        shutil.copy(outfile.name, outname)
        
        
def mmetsp_name(old_name):
    """
    Produce a new name for sequences from MMETSP
    :param old_name:
    :return:
    """
    pass

def jgi_name(old_name):
    """
    Produce a new name for sequences from JGI
    :param old_name:
    :return:
    """
    
parser = ArgumentParser(description='Unified names for diatom sequences')
parser.add_argument('-j', type=str, help='JGI FASTA')
parser.add_argument('-m', type=str, help='MMETSP FASTA')
parser.add_argument('-s', type=str, help='Synedra acus FASTA')
parser.add_argument('-i', action='store_true', help='Rename in-place')
args = parser.parse_args()

if args.j:
    rename_all(args.j, jgi_name, args.i)
if args.m:
    rename_all(args.m, mmetsp_name, args.i)
if args.s:
    rename_all(args.s, lambda x: 'SynedraAcus_'+x, args.i)
