#! /usr/bin/env python3

#  Unit tests for phylomics.py
from phylomics import BlastHit, parse_blast_line, parse_blast_file

def test_blast_lines():
    # Parsing valid line
    assert parse_blast_line('Achnanthes_exigua_ABQ43357.1	NODE_3018_length_2386_cov_4738.01	76.37	182	43	0	1	182	411	956	5e-88	 277') ==\
        BlastHit(query_id='Achnanthes_exigua_ABQ43357.1',
                 hit_id='NODE_3018_length_2386_cov_4738.01',
                 query_pos=(1, 182),
                 hit_pos=(411, 956),
                 evalue=float('5e-88'))

