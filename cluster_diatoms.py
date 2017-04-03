#! /usr/bin/env python3

from argparse import ArgumentParser
from phylome.blast_parser import *
from math import log10


def get_cluster_score(cluster, query):
    """
    Return a cluster score
    Score is a sum of hit scores, which, in turn, are -log10(evalue of best hsp)
    for every hit in query that is present in a cluster 
    :param cluster: a list of IDs
    :param query: a list of hits
    :return: 
    """
    return sum((-log10(min((x.evalue for x in hit.hsps)))
                       for hit in query if hit.hit_id in cluster))


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
        score = get_cluster_score(cluster_list[index], query)
        if score > best_score:
            best_cluster_index = index
    return best_cluster_index


parser = ArgumentParser(description='Cluster diatom sequences based on BLAST')
parser.add_argument('-b', type=str, help='BLAST tabular file')
parser.add_argument('-f', type=str, help='FASTA file')
parser.add_argument('-n', type=str, help='Base name for output files',
                    default='clusters')
args = parser.parse_args()

clusters = []
for query in iterate_by_query(parse_blast_file_to_hits(args.b)):
    best_cluster_index = get_best_cluster_index(clusters, query)
    if best_cluster_index:
        clusters[best_cluster_index].append(query[0].query_id)
    else:
        # Taking query[0] because they all have the same query ID anyway
        clusters.append([query[0].query_id])

with open('{}.clusters.list'.format(args.n), mode='w') as cluster_file:
    for index in range(len(clusters)):
        # Print cluster ID followed by the tab-separated list of sequences
        print('{}\t'.format(args.n) + '\t'.join(clusters[index]),
              file=cluster_file)
