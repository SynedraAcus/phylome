#! /usr/bin/env python3

# process outers FASTA(s): for each supplied FASTA file, iterate over records,
# removing any diatoms (ncbi taxid 2836) and checking whether there are red or
# green (NCBI taxid 2763 and 33090, respectively). Archaeplastida are written to
# {original_fasta}.archaeplastida, others are written to {original_fasta}.rest
# The script also writes a log file as the TSV. The format is:
# original_fasta    #diat   #reds   #greens #rest
#
# Also some spam in sys.stderr

from phylome.taxonomy import get_supertaxon_from_list
from argparse import ArgumentParser
import multiprocessing
import mysql.connector
import sys

def process_fasta(filename):
    """
    Process a single FASTA file. Returns a log line
    :param filename: 
    :return: 
    """
    pass

parser = ArgumentParser('Process outers list')
parser.add_argument('-f', nargs='*', type=str,
                    help='FASTA file(s)')
parser.add_argument('-l', type=str, help='log file')
parser.add_argument('-h', type=str, help='MySQL host')
parser.add_argument('-u', type=str, help='MySQL username')
parser.add_argument('-p', type=str, help='MySQL password')
args = parser.parse_args()
