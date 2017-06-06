#! /usr/bin/env python3

"""
This script clusters the sequences using seed species. While not particularly
great, this method at least definitely works. The idea is to do the following:
1. Group all the homologous sequences within the query species (whether 
inparalogs, outparalogs, whatever).
2. For each group, select all the diatom sequences that are hits to it above the
e-value cutoff.
NOTE: due to BLAST being all-against-all, step 1 happens automatically
3. For each of *those* select reference sequences. Either top-100 or all above
the e-value cutoff, whichever set is smaller.
4. Reduce the set of reference sequences (probably megataxa-optimised reduction)
5. Prepare the results in a dir or something.
"""

import sys
from argparse import ArgumentParser
from Bio import SeqIO
from phylome.blast_parser import parse_blast_file_to_hits

parser = ArgumentParser(description='Query species-based clustering')
parser.add_argument('-f', type=str, help='Query species FASTA')
parser.add_argument('-o', type=str, help='Output file')
parser.add_argument('-dia_blast', type=str, help='Diatom BLAST TSV')
parser.add_argument('-nr_blast', type=str, help='nr BLAST TSV')
parser.add_argument('-dia_evalue', type=float,
                    help='Diatom BLAST e-value cutoff')
parser.add_argument('-nr_evalue', type=float,
                    help='nr e-value cutoff')
args = parser.parse_args()

diatom_components = {}
for record in SeqIO.parse(args.f, 'fasta'):
    diatom_components[record.id] = set([record.id])

# for blast_hit in parse_blast_file_to_hits(args.dia_blast):
#     if min(x.evalue for x in blast_hit.hsps) <= args.dia_evalue:
#         for query, cluster in zip(queries, diatom_components):
#             if query == blast_hit.query_id:
#                 cluster.add(blast_hit.hit_id)
#                 break
for blast_line in open(args.dia_blast):
    arr = blast_line.split('\t')
    if float(arr[10]) <= args.dia_evalue:
        if arr[0] in diatom_components.keys():
            diatom_components[arr[0]].add(arr[1])
print('Diatom clusters built', file=sys.stderr)
with open('{}.diatom'.format(args.o), mode='w') as tmp:
    print(diatom_components, file=tmp)

nr_components = {x: set() for x in diatom_components.keys()}

for blast_line in open(args.nr_blast):
    arr = blast_line.split('\t')
    for valid_cluster in (x for x in nr_components.keys\
                          if arr[0] in diatom_components[x]):
        nr_components[valid_cluster].add(arr[1])
print('NR clusters built', file=sys.stderr)
with open('{}.nr'.format(args.o), mode='w') as tmp:
    print(nr_components, file=tmp)
# nr_components = [set() for x in diatom_components]
# for blast_hit in parse_blast_file_to_hits(args.nr_blast):
#     for index, diatom_cluster in enumerate(diatom_components):
#         if blast_hit.query_id in diatom_cluster and \
#                 min(x.evalue for x in blast_hit.hsps) <= args.nr_evalue:
#             nr_components[index].add(blast_hit.hit_id)
#
# clusters = [x.union(y) for x, y in zip(diatom_components, nr_components)]
# #  Here for debugging purposes, later here'll be a proper reducer
# with open(args.o, mode='w') as outfile:
#     for cluster in clusters:
#         print(cluster, file=outfile)
