#! /usr/bin/env python3

from argparse import ArgumentParser
from collections import Counter
from sys import stderr

from phylome.blast_parser import *


def slc(hits_graph, write_log=True, log_frequency=10000):
    """
    Return clusters in graph as detected by single-linkage clustering
    :param hits_graph: 
    :return: 
    """
    clusters = []
    cluster_id = 0
    while len(hits_graph) > 0:
        cluster = set()
        queries = [next(iter(hits_graph.keys()))]
        while len(queries) > 0:
            query = queries.pop()
            new = set(x for x in hits_graph[query] if x not in cluster and x in hits_graph)
            queries.extend(new)
            cluster.add(query)
        clusters.append(cluster)
        for x in cluster:
            if x in hits_graph:
                del hits_graph[x]
        cluster_id += 1
        if write_log and cluster_id % log_frequency == 0:
            print('Built {} clusters'.format(cluster_id), file=stderr)
    return clusters


def mcl(hits_graph, expand_factor=2, inflate_factor=1.5,
        max_loop=10, mult_factor=1):
    """
    Return clusters in hits graph as detected by MCL clustering algorithm
    Basically a wrapper around python_mcl's mcl function that allows to use
    hits graph instead of a matrix
    :param hits_graph: list of hit sets
    :return: 
    """
    # Importing here to allow slc to work without numpy
    import numpy as np
    from phylome.mcl_clustering import mcl
    keys = tuple(hits_graph.keys())
    l = len(keys)
    matrix = np.zeros((l, l))
    for x in range(len(keys)):  # IMO a bit cleaner than enumerate
        edge_probability = 1/len(hits_graph[keys[x]])
        for y in range(len(keys)):
            if keys[y] in hits_graph[keys[x]]:
                matrix[x, y] = edge_probability
    cluster_dict = mcl(matrix, expand_factor=expand_factor,
                      inflate_factor=inflate_factor, max_loop=max_loop,
                      mult_factor=mult_factor)
    del matrix  # Free some memory
    cluster_list = []
    for cluster in slc(cluster_dict):
        cluster_list.append(set(keys[i] for i in cluster))
    return cluster_list


parser = ArgumentParser(description='Cluster diatom sequences based on BLAST')
parser.add_argument('-b', type=str, help='BLAST tabular file')
parser.add_argument('-e', type=float, help='e-value cutoff. Default 1e-5',
                    default=1e-5)
parser.add_argument('-n', type=str, help='Base name for output files',
                    default='clusters')
parser.add_argument('-v', action='store_true', help='Produce STDERR output')
parser.add_argument('-c', type=int, help='Debug lines frequency. Works only with slc clustering',
                    default=10000)
parser.add_argument('-r', action='store_false', help='Ignore non-reciprocal hits. Default False.')
parser.add_argument('-a', type=str, help='Clustering algorithm. One of slc (for single-linkage clustering) or mcl (for Markov clustering).',
                    default='slc')
args = parser.parse_args()


#
if args.a not in ('slc', 'mcl'):
    stderr.write('Incorrect algorithm\n')
    quit()
# Loading data
if args.v:
    stderr.write('Loading hits... ')
hits = {}
for query in iterate_by_query(parse_blast_file_to_hits(args.b)):
    hits[query[0].query_id] = set(hit.hit_id for hit in query
                          if min(hsp.evalue for hsp in hit.hsps) < args.e)
if args.v:
    stderr.write('{} loaded.\n'.format(len(hits)))
    stderr.flush()
if args.r:
    if args.v:
        stderr.write('Removing non-reciprocal hits...')
    for query in hits.keys():
        hits[query] = set(hit for hit in hits[query]
                          if hit in hits.keys() and query in hits[hit])
    if args.v:
        stderr.write('{} remained\n'.format(len(hits)))

# Clustering
if args.v:
    stderr.write('Begin clustering\n')
if args.a == 'slc':
    clusters = slc(hits, write_log=args.v, log_frequency=args.c)
elif args.a == 'mcl':
    clusters = mcl(hits)


# Writing clusters
with open('{}.clusters.list'.format(args.n), mode='w') as cluster_file:
    for a in range(len(clusters)):
        print(a, clusters[a], file=cluster_file)

#  Writing cluster data
if args.v:
    print('Calculating cluster lengths distribution', file=stderr)
with open('{}.cluster_sizes'.format(args.n), mode='w') as data_file:
    lengths = [len(x) for x in clusters]
    lengths.sort()
    len_counts = Counter(lengths)
    for i in sorted(len_counts.keys()):
        print('{}\t{}'.format(i, len_counts[i]), file=data_file)
