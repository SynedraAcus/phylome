#! /usr/bin/env python3

from argparse import ArgumentParser

import dendropy

parser = ArgumentParser('Find red-supporting or green-supporting trees')
parser.add_argument('-n', type=str, help='Cluster names file')
parser.add_argument('-k', type=float, help='Taxa represenation', default=0.05)
args = parser.parse_args()


def paint_tree(tree, taxfile):
    """
    For each node, note how many diatoms, reds or greens are available
    :param tree:
    :param taxfile:
    :return:
    """
    global groups # Saving like a zillionth of a second on lookup
    leaf_groups = {}
    for line in open(taxfile):
        l = line.rstrip().split('\t')
        leaf_groups[l[0]] = l[1]
    # Painting leaves
    for node in tree.leaf_node_iter():
        for group in groups:
            node.annotations[group] = 0
        node.annotations[leaf_groups[node.taxon.label]] = 1
    # Painting internal nodes
    for node in tree.postorder_node_iter():
        if not node.child_nodes():
            # Skipping leaves
            continue
        for group in groups:
            node.annotations[group] = 0
        for child in node.child_node_iter():
            for group in groups:
                node.annotations[group].value += child.annotations[group].value
    
def get_tree_class(tree, k):
    """
    Process a tree, returning either 'red', 'green', 'rest' or 'neither'
    See journal, 7/XII/2017, for analysis description
    :param tree: 
    :param connection: 
    :return: 
    """
    # Tree stats
    global groups
    # Root node stats are tree stats
    counts = {x: tree.seed_node.annotations[x].value for x in groups}
    # Maximum available in alien clades
    maxes = {x: round(counts[x]*k) for x in counts}
    # Minimum necessary for a clade
    mins = {x: counts[x] - round(counts[x]*k) for x in counts}
    for node in tree.preorder_node_iter():
        if node.annotations['red'].value >= mins['red'] and \
                node.annotations['diatom'].value >= mins['diatom'] and \
                node.annotations['green'].value <= maxes['green'] and \
                node.annotations['rest'].value <= maxes['rest']:
            return 'red'
        elif node.annotations['green'].value >= mins['green'] and \
                node.annotations['diatom'].value >= mins['diatom'] and \
                node.annotations['red'].value <= maxes['red'] and \
                node.annotations['rest'].value <= maxes['rest']:
            return 'green'
        elif node.annotations['rest'].value >= mins['rest'] and \
                node.annotations['diatom'].value >= mins['diatom'] and \
                node.annotations['green'].value <= maxes['green'] and \
                node.annotations['red'].value <= maxes['red']:
            return 'rest'
    else:
        return 'neither'
    

tree_mask = 'trees/{}.dataset.fasta.dist.tsv.tre'
taxdata_mask = 'taxdata/{}.taxdata.tsv'
groups = ('red', 'green', 'rest', 'diatom')
for line in open(args.n):
    cluster_id = line.rstrip()
    tree = dendropy.Tree.get(path=tree_mask.format(cluster_id), schema='newick')
    paint_tree(tree, taxdata_mask.format(cluster_id))
    print(get_tree_class(tree, args.k))
