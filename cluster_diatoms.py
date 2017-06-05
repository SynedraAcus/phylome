#! /usr/bin/env python3

from argparse import ArgumentParser
from collections import Counter
from copy import deepcopy
from sys import stderr

from phylome.blast_parser import *


def reduce_graph(hits_graph, l):
    """
    Retain only edges between vertices that share no less than l neighbours.
    l is a proportion of shared neighbours to neighbour count of a node that has
    less neighbours
    :param hits_graph: dict of sets
    :param l: float
    :return:
    """
    # A copy of graph avoids creating artifacts by the order of edge processing
    new_graph = deepcopy(hits_graph)
    for item in hits_graph.keys():
        neighbours = hits_graph[item]
        for neighbour in neighbours:
            # Edge between the two nodes in question is considered a neighbour
            # for elimination purposes
            shared_neighbours = \
                len(hits_graph[item].intersection(hits_graph[neighbour])) + 1
            if shared_neighbours < min(len(hits_graph[item]),
                                        len(hits_graph[neighbour])) * l:
                #  Removing edges if there are not enough shared neighbours
                #  These try/except constructions are necessary in a
                #  non-reciprocal hits graph
                try:
                    new_graph[item].remove(neighbour)
                except KeyError:
                    pass
                try:
                    new_graph[neighbour].remove(item)
                except KeyError:
                    pass
    return new_graph
    

def slc(hits_graph, write_log=True, log_frequency=10000):
    """
    Return connected components in the graph
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
            print('Built {} diatom_components'.format(cluster_id), file=stderr)
    return clusters


def mcl(hits_graph, expand_factor=2, inflate_factor=1.5,
        max_loop=10, mult_factor=1):
    """
    Return diatom_components in hits graph as detected by MCL clustering algorithm
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


def red(hits_graph, l):
    """
    Reduce the graph and produce diatom_components in the reduced one
    :param hits_graph:
    :param l:
    :return:
    """
    return slc(reduce_graph(hits_graph, l))


parser = ArgumentParser(description='Cluster diatom sequences based on BLAST')
parser.add_argument('-b', type=str, help='BLAST tabular file')
parser.add_argument('-e', type=float, help='e-value cutoff. Default 1e-5',
                    default=1e-5)
parser.add_argument('-n', type=str, help='Base name for output files',
                    default='diatom_components')
parser.add_argument('-v', action='store_true', help='Produce STDERR output')
parser.add_argument('-c', type=int, help='Debug lines frequency. Works only with slc clustering',
                    default=10000)
parser.add_argument('-r', action='store_false', help='Ignore non-reciprocal hits. Default False.')
parser.add_argument('-a', type=str, help='Clustering algorithm. One of slc (for single-linkage clustering) or mcl (for Markov clustering) or red fpr graph reduction.',
                    default='slc')
parser.add_argument('-i', type=float, default=1.5,
                    help='Inflation factor for MCL algorithm. Default 1.5')
parser.add_argument('-l', type=float, default=0.5,
                    help='Shared neighbours proportion for graph reduction')
args = parser.parse_args()


#
if args.a not in ('slc', 'mcl', 'red'):
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
    clusters = mcl(hits, inflate_factor=args.i)
elif args.a == 'red':
    clusters = red(hits, l=args.l)



# Writing diatom_components
with open('{}.diatom_components.list'.format(args.n), mode='w') as cluster_file:
    for a in range(len(clusters)):
        print(a, clusters[a], file=cluster_file)

#  Writing diatom_cluster data
if args.v:
    print('Calculating diatom_cluster lengths distribution', file=stderr)
with open('{}.cluster_sizes'.format(args.n), mode='w') as data_file:
    lengths = [len(x) for x in clusters]
    lengths.sort()
    len_counts = Counter(lengths)
    for i in sorted(len_counts.keys()):
        print('{}\t{}'.format(i, len_counts[i]), file=data_file)
