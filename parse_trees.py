#! /usr/bin/env python3

import dendropy
from argparse import ArgumentParser
parser = ArgumentParser('Find red-supporting or green-supporting trees')
parser.add_argument('-t', type=str, nargs='*', help='Newick trees')
parser.add_argument('-o', type=str, help='MySQL host')
parser.add_argument('-d', type=str, help='MySQL database')
parser.add_argument('-u', type=str, help='MySQL username')
parser.add_argument('-p', type=str, help='MySQL password')
args = parser.parse_args()

def get_tree_class(tree, cursor):
    """
    Process a tree, returning either 'red', 'green' or 'neither'
    See journal, 7/XII/2017, for analysis description
    :param tree: 
    :param connection: 
    :return: 
    """
    pass

for treefile in args.t:
    tree = dendropy.Tree.get(path=treefile, schema='newick')
    tree_class = get_tree_class(tree, cursor)

