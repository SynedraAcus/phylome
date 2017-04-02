#! /usr/bin/env python3

from argparse import ArgumentParser

from Bio import SeqIO


def rename_sequence(record, renamer=lambda x: x):
    # record.description = renamer(record.description)
    record.id = renamer(record.description)
    record.description = ''
    record.name = ''
    #  return for readability only
    return record


def rename_all(file, renamer):
    # Iterators are the best thing since sliced bread
    SeqIO.write((rename_sequence(x, renamer) for x in SeqIO.parse(file, format='fasta')),
                handle='{}.renamed'.format(file),
                format='fasta')
        
        
def mmetsp_name(old_name):
    """
    Produce a new name for sequences from MMETSP
    :param old_name:
    :return:
    """
    global mmetsp_regex
    l = old_name.split(' ')
    name = l[0]
    for x in l:
        if 'ORGANISM' in x:
            l2 = old_name.split('\"')
            name = l2[1] + '_' + name


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
args = parser.parse_args()
mmetsp_regex = re.compile('/ORGANISM="(\w+)"')
if args.j:
    rename_all(args.j, jgi_name, args.i)
if args.m:
    rename_all(args.m, mmetsp_name, args.i)
if args.s:
    rename_all(args.s, renamer=lambda x: 'Synedra_Acus_'+x)
