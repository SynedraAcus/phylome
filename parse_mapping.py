#! /usr/bin/env python3

# p

from argparse import ArgumentParser
import re

parser = ArgumentParser('Process the IQ-tree likelihood mapping data')
parser.add_argument('-f', nargs='*', help='iqtree files')
args = parser.parse_args()

clusters = ('diatom', 'red', 'green', 'rest')
cutoff_line='-'*144+'\n'
regex = re.compile('\(([\d.]+) *\)')
segment_values = []
for iqfile in args.f:
    # Read the summary segment support from seven-segment table
    # It'd be right after the table, which begins and ends with cutoff line
    wait = 0
    for line in open(iqfile):
        if wait == 2:
            segment_values.append([float(x) for x in re.findall(regex, line)])
            break
        if line == cutoff_line:
            wait += 1

summary = [0 for x in range(7)]
for cluster_run in segment_values:
    for index, value in enumerate(cluster_run):
        summary[index] += value
print(summary)
