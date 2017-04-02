#! /usr/bin/env python3

import shutil
import sys
from argparse import ArgumentParser
from tempfile import NamedTemporaryFile

from Bio import SeqIO

parser = ArgumentParser(description='Filter aminoacid sequences in FASTA file by length and percentage of X\'es')
parser.add_argument('-f', type=str, nargs='+', help='FASTA file(s)')
parser.add_argument('-x', type=float, default=0.1,
                    help='Maximum allowed percentage of X\'s')
parser.add_argument('-l', type=int, help='Minimum allowed length. Default 100',
                    default=100)
parser.add_argument('-i', action='store_true',
                    help='Filter in-place. Disabled by default')
parser.add_argument('-v', action='store_true',
                    help='Verbose output')
args = parser.parse_args()


for fasta_file in args.f:
    with NamedTemporaryFile(mode='w+') as output_handle:
        record_count = 0
        accepted_record_count = 0
        for record in SeqIO.parse(fasta_file, format='fasta'):
            record_count += 1
            if len(record) > args.l:
                x_perc = (record.seq.count('X') + record.seq.count('x'))/len(record)
                if x_perc < args.x:
                    accepted_record_count += 1
                    SeqIO.write(record, handle=output_handle, format='fasta')
            if args.v and record_count % 10000 == 0:
                sys.stderr.write('Accepted {0} out of {1} records ({2}%) from {3}\n'.\
                    format(accepted_record_count, record_count,
                    int(100*accepted_record_count/record_count), fasta_file))
        if args.v:
            sys.stderr.write('Accepted {0} out of {1} records ({2}%) from {3}\n'.format(
                accepted_record_count, record_count,
                int(100*accepted_record_count/record_count), fasta_file
            ))
            sys.stderr.flush()
        if args.i:
            shutil.copy(output_handle.name, fasta_file)
        else:
            dest = '{}.filtered'.format(fasta_file)
            shutil.copy(output_handle.name, dest)
