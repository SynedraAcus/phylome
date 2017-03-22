#! /usr/bin/env python3

#  Unit tests for blast_parser.py
import pytest
from blast_parser import BlastHSP, parse_blast_line, parse_blast_file_to_hsps, \
    assemble_hits, BlastHit, parse_blast_file_to_hits


def test_valid_blast_lines():
    # Parsing valid line
    assert parse_blast_line('Achnanthes_exigua_ABQ43357.1	NODE_3018_length_2386_cov_4738.01	76.37	182	43	0	1	182	411	956	5e-88	 277') == \
           BlastHSP(query_id='Achnanthes_exigua_ABQ43357.1',
                    hit_id='NODE_3018_length_2386_cov_4738.01',
                    query_pos=(1, 182),
                    hit_pos=(411, 956),
                    evalue=float('5e-88'))
    # Numeric IDs?
    assert parse_blast_line('1	2	76.37	182	43	0	1	182	411	956	5e-88	 277') == \
           BlastHSP(query_id='1',
                    hit_id='2',
                    query_pos=(1, 182),
                    hit_pos=(411, 956),
                    evalue=float('5e-88'))


def test_invalid_blast_lines():
    # A single field missing
    with pytest.raises(ValueError):
        parse_blast_line('1	2	76.37	182	43	0	1	411	956	5e-88	 277')
    # One field too much
    with pytest.raises(ValueError):
        parse_blast_line('0\t1	2	76.37	182	43	0	1	182	411	956	5e-88	 277')
    with pytest.raises(ValueError):
        parse_blast_line('')
    with pytest.raises(ValueError):
        parse_blast_line('Foo bar baz')

# Calling list() everywhere because generator is not evaluated until someone
# actually needs a value


def test_valid_blast_file_to_hsps():
    #  Check that the generator returns 50 BlastHSP`s without raising anything
    # A filename passed as a string
    l = list(parse_blast_file_to_hsps('test_data/BLAST_test.tsv'))
    assert len(l) == 50 and all(isinstance(x, BlastHSP) for x in l)
    # A filename passed as a filehandle
    l = list(parse_blast_file_to_hsps(open('test_data/BLAST_test.tsv')))
    assert len(l) == 50 and all(isinstance(x, BlastHSP) for x in l)


def test_invalid_blast_file_to_hsps():
    with pytest.raises(TypeError):
        g = list(parse_blast_file_to_hsps(42))
    with pytest.raises(TypeError):
        g = list(parse_blast_file_to_hsps([1, 2, 'qwerty']))
    # File doesn't exist:
    with pytest.raises(FileNotFoundError):
        g = list(parse_blast_file_to_hsps('asdf'))
    # A file exists, but isn't BLAST
    with pytest.raises(ValueError):
        g = list(parse_blast_file_to_hsps('README.md'))
    # Same for handle
    with pytest.raises(ValueError):
        g = list(parse_blast_file_to_hsps(open('README.md', mode='r')))

#  Similar to previous two


def test_valid_blast_file_to_hits():
    #  Check that the generator returns 45 BlastHit`s without raising anything
    # A filename passed as a string
    l = list(parse_blast_file_to_hits('test_data/BLAST_test.tsv'))
    assert len(l) == 45 and all(isinstance(x, BlastHit) for x in l)
    # A filename passed as a filehandle
    l = list(parse_blast_file_to_hits(open('test_data/BLAST_test.tsv')))
    assert len(l) == 45 and all(isinstance(x, BlastHit) for x in l)


def test_invalid_blast_file_to_hits():
    with pytest.raises(TypeError):
        g = list(parse_blast_file_to_hits(42))
    with pytest.raises(TypeError):
        g = list(parse_blast_file_to_hits([1, 2, 'qwerty']))
    # File doesn't exist:
    with pytest.raises(FileNotFoundError):
        g = list(parse_blast_file_to_hits('asdf'))
    # A file exists, but isn't BLAST
    with pytest.raises(ValueError):
        g = list(parse_blast_file_to_hits('README.md'))
    # Same for handle
    with pytest.raises(ValueError):
        g = list(parse_blast_file_to_hits(open('README.md', mode='r')))


def test_assemble_hits():
    # This file contains 45 hits, only 5 out of which have 2 hsps and the rest 1
    hits = list(assemble_hits(parse_blast_file_to_hsps('test_data/BLAST_test.tsv')))
    assert len(hits) == 45
    assert all(isinstance(x, BlastHit) for x in hits)
    assert len([x for x in hits if len(x.hsps) == 2]) == 5
    assert len([x for x in hits if len(x.hsps) > 2 or len(x.hsps) < 1]) == 0

