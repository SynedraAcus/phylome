#! /usr/bin/env python3

from argparse import ArgumentParser
from math import log10
from sys import stderr

from phylome.blast_parser import *


def get_hit_score(hit):
    """
    Score is -log10(evalue) or 1000 if evalue==0
    :param evalue: 
    :return: 
    """
    return max(-log10(hsp.evalue) if hsp.evalue > 0 else 1000 for hsp in hit.hsps)


def get_best_cluster_index(cluster_list, query):
    """
    Return the index of the best cluster or None if no cluster is present
    :param cluster_list: 
    :param query: 
    :return: 
    """
    best_cluster_index = None
    best_score = 0
    for index in range(len(cluster_list)):
        # Iterate once instead of checking for the presence of *every* hit on
        # the next line. Maybe will save some time with large queries being
        # checked against the wrong clusters
        if not any(x.hit_id in cluster_list[index] for x in query):
            continue
        score = sum(get_hit_score(hit) for hit in query
                    if hit.hit_id in cluster_list[index])
        if score > best_score:
            best_cluster_index = index
    return best_cluster_index


parser = ArgumentParser(description='Cluster diatom sequences based on BLAST')
parser.add_argument('-b', type=str, help='BLAST tabular file')
parser.add_argument('-f', type=str, help='FASTA file')
parser.add_argument('-n', type=str, help='Base name for output files',
                    default='clusters')
parser.add_argument('-v', action='store_true', help='Produce STDERR output')
parser.add_argument('-c', type=int, help='Debug lines frequency')
parser.add_argument('-t', action='store_true', help='Store clusters file every 10k queries')
args = parser.parse_args()

# Loading data
hits = {}
for query in iterate_by_query(parse_blast_file_to_hits(args.b)):
    hits[query[0].query_id] = set(hit.hit_id for hit in query)
print('Removing non-reciprocal hits', file=stderr)
stderr.flush()
for query in hits.keys():
    hits[query] = set(hit for hit in hits[query] if query in hits[hit])

# Clustering
clusters = []
cluster_id = 0
with open('{}.clusters.list'.format(args.n), mode='w') as cluster_file:
    while len(hits) > 0:
        cluster = set()
        queries = [next(iter(hits.keys()))]
        while len(queries) > 0:
            query = queries.pop()
            new = set(x for x in hits[query] if x not in cluster)
            queries.extend(new)
            cluster.add(query)
        clusters.append(cluster)
        cluster_id += 1
        print(cluster_id, cluster, file=cluster_file)
        for x in cluster:
            del hits[x]
