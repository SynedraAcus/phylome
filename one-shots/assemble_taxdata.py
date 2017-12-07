#! /usr/bin/env python3

from argparse import ArgumentParser
from Bio import SeqIO
import os


def get_path(mask, cluster_id):
    """
    Return mask.format(cluster_id)+'.reduced', or, if that doesn't exist,
    mask.format(cluster_id), or, failing that, None
    :param mask: 
    :param cluster_id: 
    :return: 
    """
    if os.path.isfile(mask.format(cluster_id)+'.reduced'):
        return mask.format(cluster_id)+'.reduced'
    elif os.path.isfile(mask.format(cluster_id)):
        return mask.format(cluster_id)
    else:
        return None


parser = ArgumentParser('Assemble taxonomic mappings')
parser.add_argument('-l', type=str, help='List file')
parser.add_argument('-d', type=str, help='Output directory')
args = parser.parse_args()

diatom_mask = '/home/amorozov/phylome_data/clustering/groups/{}.diatoms.fasta'
# These should be checked first at {}.reduced, then, if that failed, at just
# mask
red_mask = '/home/amorozov/phylome_data/clustering/outers/{}.external.fasta.2763'
green_mask = '/home/amorozov/phylome_data/clustering/groups/{}.external.fasta.33090'
rest_mask = '/home/amorozov/phylome_data/clustering/groups/{}.external.fasta.rest'
if not os.path.isdir(args.d):
    os.mkdir(args.d)
for line in args.l:
    cluster_id = line.rstrip()
    with open(os.path.join(args.d, '{}.taxdata.tsv'.format(cluster_id)),
              mode='w') as outhandle:
        for record in SeqIO.parse(diatom_mask.format(cluster_id), 'fasta'):
            print('{}\tdiatom'.format(record.id), file=outhandle)
        # Checking masks
        # Could've made that a function, but whatever
        red_path = get_path(red_mask, cluster_id)
        if red_path:
            for record in SeqIO.parse(red_path, 'fasta'):
                print('{}\tred'.format(record.id), file=outhandle)
        green_path = get_path(green_mask)
        if green_path:
            for record in SeqIO.parse(green_path, 'fasta'):
                print('{}\tgreen'.format(record.id), file=outhandle)
        rest_path = get_path(rest_mask)
        if rest_path:
            for record in SeqIO.parse(rest_path, 'fasta'):
                print('{}\trest'.format(record.id), file=outhandle)
