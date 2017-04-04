#! /usr/bin/env python3

from argparse import ArgumentParser
from phylome.blast_parser import *
from math import log10
from sys import stderr

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
parser.add_argument('-v', type=str, help='Produce STDERR output')
parser.add_argument('-c', type=int, help='Debug lines frequency')
args = parser.parse_args()

clusters = []
query_count = 0
for query in iterate_by_query(parse_blast_file_to_hits(args.b)):
    best_cluster_index = get_best_cluster_index(clusters, query)
    if best_cluster_index:
        clusters[best_cluster_index].append(query[0].query_id)
    else:
        # Taking query[0] because they all have the same query ID anyway
        clusters.append([query[0].query_id])
    query_count += 1
    if query_count % args.c == 0:
        stderr.write('{} queries processed, {} clusters created'.format(
            query_count, len(clusters)
        ))

with open('{}.clusters.list'.format(args.n), mode='w') as cluster_file:
    for index in range(len(clusters)):
        # Print cluster ID followed by the tab-separated list of sequences
        print('{}\t'.format(str(index)) + '\t'.join(clusters[index]),
              file=cluster_file)
