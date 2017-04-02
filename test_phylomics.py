#! /usr/bin/env python3

# Not my code, but needed for tests
import mysql.connector
import pytest

# Stuff to be tested
from phylome.blast_parser import BlastHSP, parse_blast_line, \
    parse_blast_file_to_hsps, \
    assemble_hits, BlastHit, parse_blast_file_to_hits, iterate_by_query
from phylome.multiplicates import is_duplicate
from phylome.taxonomy import descend_taxon_tree, get_taxa_list, is_taxon_member, \
    get_supertaxon_from_list


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
    assert len(l) == 46 and all(isinstance(x, BlastHit) for x in l)
    # A filename passed as a filehandle
    l = list(parse_blast_file_to_hits(open('test_data/BLAST_test.tsv')))
    assert len(l) == 46 and all(isinstance(x, BlastHit) for x in l)


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
    # This file contains 46 hits, only 5 out of which have 2 hsps and the rest 1
    hits = list(assemble_hits(parse_blast_file_to_hsps('test_data/BLAST_test.tsv')))
    assert len(hits) == 46
    assert all(isinstance(x, BlastHit) for x in hits)
    assert len([x for x in hits if len(x.hsps) == 2]) == 4
    assert len([x for x in hits if len(x.hsps) > 2 or len(x.hsps) < 1]) == 0
    

def test_iterate_by_query():
    # The same file as in previous test. It contains 9 unique queries
    l = list(iterate_by_query(parse_blast_file_to_hits('test_data/BLAST_test.tsv')))
    assert len(l) == 9
    # These queries (in order that they appear) have this many hits
    assert [len(x) for x in l] == [5, 3, 5, 6, 6, 6, 6, 4, 5]


def test_duplicates():
    # Test the duplicate-detection methods.
    # query1 is duplicate, query2 isn't, single_hsp is guess what
    duplicate_hit, non_duplicate_hit, single_hsp = assemble_hits(
        [BlastHSP(query_id='query1', hit_id='hit',
                  evalue=0.0001, query_pos=(1, 100),
                  hit_pos=(1, 100)),
         BlastHSP(query_id='query1', hit_id='hit',
                  evalue=0.0001, query_pos=(150, 250),
                  hit_pos=(5, 95)),
         BlastHSP(query_id='query2', hit_id='hit',
                  evalue=0.0001, query_pos=(1, 100),
                  hit_pos=(1, 100)),
         BlastHSP(query_id='query2', hit_id='hit',
                  evalue=0.001, query_pos=(120, 170),
                  hit_pos=(110, 180)),
         BlastHSP(query_id='query3', hit_id='hit',
                  evalue=0.0001, query_pos=(1, 100),
                  hit_pos=(1, 100))])
    assert not is_duplicate(single_hsp)
    assert is_duplicate(duplicate_hit)
    assert not is_duplicate(non_duplicate_hit)


def test_taxa_descent():
    #  All mySQL tests are kept in a single function to avoid lengthy creation
    #  of a connection for every test I need to run.
    #  This test will only run if you have bioSQL mySQL DB running with default
    #  username and password. You probably shouldn't.
    connection = mysql.connector.connect(host='localhost', user='root',
                                         password='password', database='biosql')
    cursor = connection.cursor()
    #  An ancestor chain for a valid taxon, in this case genus synedra
    # assert list(descend_taxon_tree(156135, cursor)) == [
    #     16060, 16059, 16058, 16057, 2231, 15888, 2166, 101881]
    assert list(descend_taxon_tree(191584, cursor)) == [
        33856, 33855, 33854, 33853, 2836, 33634, 2759, 131567]
    #  ValueError on the invalid taxon_id
    with pytest.raises(ValueError):
        list(descend_taxon_tree(10000000, cursor))
        list(descend_taxon_tree(0, cursor))
    with pytest.raises(TypeError):
        list(descend_taxon_tree('foo', cursor))
    #  The same in reverse, and as list, not generator
    assert get_taxa_list(191584, cursor) == [
        131567, 2759, 33634, 2836, 33853, 33854, 33855, 33856]
    #  Correct ancestor
    assert is_taxon_member(191584, 2836, cursor)
    #  Diatoms are definitely not bacteria
    assert not is_taxon_member(191584, 2, cursor)
    #  Again, are they bacteria, archaea or eukaryotes?
    assert get_supertaxon_from_list(191584, [2, 2157, 2759], cursor) == 2759
    #  Are diatoms mammals or birds?
    with pytest.raises(ValueError):
        get_supertaxon_from_list(191584, [40674, 8782], cursor)
    
