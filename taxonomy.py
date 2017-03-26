"""
A bunch of functions for making a checking taxonomy statements.

All functions in this file rely on NCBI taxonomy database being available via
mySQL and require a mysql.connector.cursor as one of the arguments. All these
things can probably be written better using pure SQL, but I don't know SQL and
I kinda need these things right now.
"""


def descend_taxon_tree(starting_taxon, cursor):
    """
    A generator that takes a starting taxon and yields its ancestors until it
    hits the root.
    :param starting_taxon:
    :param cursor:
    :return:
    """
    pass

    
def get_taxa_list(taxon_id, cursor):
    """
    Return a list of taxa for a given taxon_id.
    Accepts a taxon id and returns a list of all taxa it belongs (as NCBI taxon
    ids), sorted from larger taxa to smaller ones. Does not include root
    :param taxon_id:
    :param cursor:
    :return:
    """


def is_taxon_member(taxon1, taxon2, cursor):
    """
    Return True if taxon 1 is a subtaxon of taxon2, False otherwise
    :param taxon_id:
    :param cursor:
    :return:
    """
    pass
    

def get_supertaxon_from_list(taxon, taxa_list, cursor):
    """
    Return a taxon from a given list that is a supertaxon of a given taxon.
    Accepts a taxon id and a list of taxon ids. If one of the taxa in the list
    is a supertaxon of the query taxon, returns the id of that one. If neither
    is, ie the function has descended to the root without hitting one of the
    query taxa, raises ValueError.
    If taxa in taxa_list are nested, the one nearest to the query taxon is
    returned.
    :param taxon:
    :param taxa_list:
    :param cursor:
    :return:
    """
    pass
