#! /usr/bin/env python3

import collections.abc
from argparse import ArgumentParser
from Bio import SeqIO
from kmers.kmers import Composition, ffp_distance


class Element:
    """A simple data class: a sequence and its composition."""
    def __init__(self, sequence=None, k=3):
        self.sequence = sequence
        self.id = sequence.id
        self.composition = Composition(k=k, seq=sequence)


class Cluster(collections.abc.MutableSet):
    """
    A cluster of Sequence objects.
    """
    def __init__(self, elements, cutoff=0.3):
        assert all(isinstance(x, Element) for x in elements)
        self.sequences = set(elements)
        self.cutoff = cutoff

    def __contains__(self, item):
        if item in self.sequences:
            return True

    def __iter__(self):
        return self.sequences.__iter__()

    def __len__(self):
        return len(self.sequences)

    def should_accept(self, item):
        """
        Return True if the item is similar enough to all items in cluster
        :param item:
        :return:
        """
        for a in self:
            # No problem if it walks through this zero times
            if ffp_distance(a.composition, item.composition) > self.cutoff:
                return False
        return True

    def add(self, item):
        assert isinstance(item, Element)
        self.sequences.add(item)

    def discard(self, value):
        self.sequences.discard(value)

if __name__ == '__main__':
    parser = ArgumentParser(description='A fast k-mer based clusterization tool')
    parser.add_argument('-f', type=str, nargs='*', help='FASTA file(s)')
    parser.add_argument('-t', type=float, default=0.3,
                        help='Clustering distance cutoff')
    args = parser.parse_args()

    clusters = []
    for fasta_file in args.f:
        for record in SeqIO.parse(fasta_file, format='fasta'):
            element = Element(record)
            accepted = False
            for cluster in clusters:
                if cluster.should_accept(element):
                    cluster.add(element)
                    accepted = True
            if not accepted:
                clusters.append(Cluster((element, )))
                print('.')
    print('Built {0} clusters'.format(len(clusters)))
    for cluster in clusters:
        print(' '.join((x.id for x in cluster)))

