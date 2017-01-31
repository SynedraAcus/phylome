#! /usr/bin/env python3

import collections.abc
from argparse import ArgumentParser
from Bio import SeqIO
from kmers.kmers import Composition, ffp_distance


class Sequence():
    """A simple data class: a sequence and its composition."""
    def __init__(self, sequence=None, k=3):
        self.sequence = sequence
        self.composition = Composition(k=k, seq=sequence)


class Cluster(collections.abc.MutableSet):
    """
    A cluster of Sequence objects
    """
    def __init__(self):
        pass

    def __contains__(self, item):
        pass

    def __iter__(self):
        pass

    def __len__(self):
        pass

    def add(self, item):
        pass

    def discard(self, value):
        pass

if __name__ == 'main':
    parser = ArgumentParser(description='A fast k-mer based clusterization tool')
    parser.add_argument('-f', type=str, nargs='*', help='FASTA file(s)')
    parser.add_argument('-t', type=float, default=0.3,
                        help='Clustering distance cutoff')
    clusters = set()

