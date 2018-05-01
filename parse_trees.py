#! /usr/bin/env python3

from argparse import ArgumentParser

import dendropy

parser = ArgumentParser('Find red-supporting or green-supporting trees')
parser.add_argument('-n', type=str, help='Cluster names file')
parser.add_argument('-k', type=float, help='Taxa represenation', default=0.05)
parser.add_argument('-s', action='store_true', help='Run sister analysis')
parser.add_argument('-c', action='store_true', help='Consider sister size for sister analysis')
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
    # Minimum necessary for a clade
    mins = {x: counts[x] - round(counts[x]*k) for x in counts}
    for node in tree.postorder_node_iter():
        if node.annotations['red'].value >= mins['red'] and \
                node.annotations['diatom'].value >= mins['diatom'] and \
                node.annotations['green'].value <= mins['green'] and \
                node.annotations['rest'].value <= mins['rest']:
            return 'red'
        elif node.annotations['green'].value >= mins['green'] and \
               node.annotations['diatom'].value >= mins['diatom'] and \
               node.annotations['red'].value <= mins['red'] and \
               node.annotations['rest'].value <= mins['rest']:
            return 'green'
        elif node.annotations['rest'].value >= mins['rest'] and \
                node.annotations['diatom'].value >= mins['diatom'] and \
                node.annotations['green'].value <= mins['green'] and \
                node.annotations['red'].value <= mins['red']:
            return 'rest'
    else:
        return 'neither'
    
def sister_analysis(tree, consider_sister):
    """
    For every diatom-only clade sister to nondiatoms,
    who precisely is it sister to?
    """
    global groups
    node_counts = {x:0 for x in groups}
    totals = {x: tree.seed_node.annotations[x].value for x in groups}
    for node in tree.postorder_internal_node_iter():
        sis = None
        weight = None
        for child in node.child_node_iter():
            if child.annotations['diatom'].value == 0:
                # This is a nondiatom node
                # sis can be harmlessly overwritten if a diatom node is not encountered
                sis = {x: child.annotations[x].value for x in groups}
            elif all([child.annotations[x].value == 0 for x in ('red', 'green', 'rest')]):
                # This is a diatom node
                weight = child.annotations['diatom'].value/totals['diatom']
        if sis and weight:
            nonzeros = [x for x in sis if sis[x] > 0]
            if len(nonzeros) == 1:
                if consider_sister:
                    node_counts[nonzeros[0]] += weight*(sis[nonzeros[0]]/totals[nonzeros[0]])
                else:
                    node_counts[nonzeros[0]] += weight
    return node_counts
                

taxdata_mask = 'taxdata/{}.clean.tsv'
tree_mask = 'ml/{}.dataset.fasta.aln.fasta.treefile'
groups = ('red', 'green', 'rest', 'diatom')
running_weights = {x: 0 for x in groups}
for line in open(args.n):
    cluster_id = line.rstrip()
    tree = dendropy.Tree.get(path=tree_mask.format(cluster_id),
                schema='newick', preserve_underscores=True)
    paint_tree(tree, taxdata_mask.format(cluster_id))
    if args.s:
        c = sister_analysis(tree, args.c)
        for x in groups:
            running_weights[x] += c[x]
    else:
        print(get_tree_class(tree, args.k))
if args.s:
    print(running_weights)
