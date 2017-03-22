#! /usr/bin/env python3

"""A collection of functions and methods for phylomic analysis"""

from collections import namedtuple
from io import TextIOBase

# Lightweight BLAST table parser
BlastHSP = namedtuple('BlastHSP', ['query_id', 'hit_id', 'query_pos', 'hit_pos',
                                   'evalue'])
BlastHit = namedtuple('BlastHit', ['query_id', 'hit_id', 'hsps'])


def parse_blast_line(line):
    """
    Take a single line from BLAST file and return a BlastHSP.
    It assumes the default tabular output from BLAST 2.2.28+ (generated with
    `-outfmt 6` without any additional format specifications). May or may not
    work with other BLAST versions.
    :param line: str
    :return: BlastHSP
    """
    line = line.rstrip('\n')
    arr = line.split('\t')
    if not len(arr) == 12:
        raise ValueError('Incorrect BLAST line {}'.format(line))
    #  There is no explicit validity check, but hopefully IndexError or some
    #  exception of the type converters will be raised on the invalid string
    return BlastHSP(query_id=arr[0],
                    hit_id=arr[1],
                    query_pos=(int(arr[6]), int(arr[7])),
                    hit_pos=(int(arr[8]), int(arr[9])),
                    evalue=float(arr[10]))


def parse_blast_file_to_hsps(filename):
    """
    Take a BLAST output file and generate BlastHSP instances from its lines.
    It assumes the default tabular output from BLAST 2.2.28+ (generated with
    `outfmt 6` or `outfmt 7` without any additional format specifications).
    Any line starting with `#` is considered a comment and omitted.
    `file` can be either a filehandle in text mode or a filename. In latter case
    it would be opened in `'r'` mode.
    :param filename: str or filehandle
    :return:
    """
    if isinstance(filename, str):
        handle = open(filename, mode='r')
    elif isinstance(filename, TextIOBase) and filename.readable():
        handle = filename
    else:
        raise TypeError('Only string or readable text-mode filehandle accepted by parse_blast_file')
    for line in handle:
        if not line[0] == '#':
            yield parse_blast_line(line)
    # Not closing a filehandle because it can be used by the calling code.


#  Assemble hits
def assemble_hits(iterable):
    """
    Groups HSPs into BlastHit instances
    Takes an iterable of BlastHSPs, yields BlastHits
    :param iterable:
    :return:
    """
    y = []
    current_hit = ''
    current_query = ''
    for blast_hsp in iterable:
        #  The following two checks will only run when the generator is created
        if not current_hit:
            current_hit = blast_hsp.hit_id
        if not current_query:
            current_query = blast_hsp.query_id
        y.append(blast_hsp)
        if current_hit !=blast_hsp.hit_id or current_query != blast_hsp.query_id:
            yield BlastHit(query_id=current_query, hit_id=current_hit, hsps=y)
            y = []
            current_hit = blast_hsp.hit_id
            current_query = blast_hsp.query_id
