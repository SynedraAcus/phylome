#! /usr/bin/env python3

"""
This script clusters the sequences using seed species. While not particularly
great, this method at least definitely works. The idea is to do the following:
1. Group all the homologous sequences within the query species (whether 
inparalogs, outparalogs, whatever.
2. For each group, select all the diatom sequences that are hits to it above the
e-value cutoff.
3. For each of *those* select reference sequences. Either top-100 or all above
the e-value cutoff, whichever set is smaller.
4. Reduce the set of reference sequences (probably megataxa-optimised reduction)
5. Prepare the results in a dir or something.
"""

from argparse import ArgumentParser
from Bio import SeqIO
from phylome.blast_parser import parse_blast_file_to_hits

parser = ArgumentParser(description='Query species-based clustering')
parser.add_argument('-dia_blast', type=str, help='Diatom BLAST TSV')
parser.add_argument('-nr_blast', type=str, help='nr BLAST TSV')
parser.add_argument('-dia_evalue', type=float,
                    help='Diatom BLAST e-value cutoff')
parser.add_argument('-nr_evalue', type=float,
                    help='nr e-value cutoff')
args = parser.parse_args()

