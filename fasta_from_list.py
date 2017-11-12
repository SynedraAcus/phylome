#! /usr/bin/env python3

import os
import sys
from argparse import ArgumentParser

from Bio import SeqIO

from phylome.blast_parser import parse_blast_file_to_hits


def batches(d, batch_size):
    """
    Yields subdicts of a given size (so many key/value pairs)
    :param dict:
    :param batch_size:
    :return:
    """
    r = {}
    c = 0
    for key in d.keys():
        r[key] = d[key]
        c += 1
        if c == batch_size:
            yield r
            r = {}
            c = 0
    yield r
    return


parser = ArgumentParser(description='Generate a lot of FASTAs using list file')
parser.add_argument('-l', type=str, help='List file')
parser.add_argument('-f', type=str, help='FASTA file(s)')
parser.add_argument('-o', type=str, help='Output basename. Default `cluster`',
                    default='cluster')
parser.add_argument('-d', type=str, default='groups',
                    help='Output directory to be created')
parser.add_argument('-b', type=str,
                    help='BLAST file of clustered seqs against other database')
parser.add_argument('--db', type=str,
                    help='External sequences FASTA')
parser.add_argument('--min', type=int, default=30,
                    help='Minimal cluster size')
parser.add_argument('--max', type=int, default=200,
                    help='Maximal cluster size')
parser.add_argument('--min_species', type=int, default=10,
                    help='Minimum species count in cluster')
parser.add_argument('--evalue', type=float, default=1e-30,
                    help='Evalue cutoff for external sequences')
parser.add_argument('--batch', type=int, help='Batch size', default=1000)
args = parser.parse_args()

clusters = {}
counter = 0
os.mkdir(args.d)
# Assumes list to consist of tab-separated ID lists, one cluster in each line
print('Loading clusters from {}'.format(args.l), flush=True, file=sys.stderr)
for line in open(args.l):
    ids = line.rstrip().split('\t')
    species = len({x.split('|')[0] for x in ids})
    if args.min <= len(ids) <= args.max and species >= args.min_species:
        clusters[counter] = set(ids)
        counter += 1
print('Loaded {} clusters'.format(counter), flush=True, file=sys.stderr)

print('Loading diatom FASTAs', flush=True, file=sys.stderr)
for cluster_batch in batches(clusters, args.batch):
    # Creating output filehandles
    diatom_filehandles = {x: open('{0}/{1}{2}.diatoms.fasta'.format(args.d,
                                                                    args.o,
                                                                    x),
                                  mode='w+')
                          for x in cluster_batch.keys()}
    for record in SeqIO.parse(args.f, 'fasta'):
        for index, cluster in cluster_batch.items():
            if record.id in cluster:
                SeqIO.write(record, diatom_filehandles[index], 'fasta')
    print('Written a batch', file=sys.stderr, flush=True)

if not args.b:
    # No external sequences
    quit()
    
# Assembling BLAST hit lists
print('Parsing BLAST file {}'.format(args.b), flush=True, file=sys.stderr)
other_seqs = {x: set() for x in clusters.keys()}
for hit in parse_blast_file_to_hits(args.b):
    for cluster_id, cluster in clusters.items():
        if hit.query_id in cluster:
            other_seqs[cluster_id].add(hit.hit_id)

print('Parsing external FASTA {}'.format(args.db), flush=True, file=sys.stderr)
# Sequences absent from the FASTA, but present in BLAST, are silently ignored
for batch in batches(other_seqs, args.batch):
    external_filehandles = {x: open('{0}/{1}{2}.external.fasta'.format(args.d,
                                                                       args.o,
                                                                       x),
                                    mode='w+')
                          for x in batch.keys()}
    for record in SeqIO.parse(args.db, 'fasta'):
        for cluster_id, cluster in batch.items():
            if record.id in cluster:
                SeqIO.write(record, external_filehandles[cluster_id])
print('Done', flush=True, file=sys.stderr)
