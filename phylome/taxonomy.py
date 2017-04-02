"""
A bunch of functions for making and checking taxonomy statements.

All functions in this file rely on NCBI taxonomy database being available via
mySQL and require a mysql.connector.cursor as one of the arguments. All these
things can probably be written trivially using pure SQL, but I don't know SQL
and I kinda need them right now.
"""


def descend_taxon_tree(starting_taxon, cursor):
    """
    A generator that takes a starting taxon and yields its ancestors until it
    hits the root.
    Returns neither query taxon nor root. Raises ValueError if supplied an
    invalid taxon ID.
    :param starting_taxon:
    :param cursor:
    :return:
    """
    if not isinstance(starting_taxon, int):
        raise TypeError('Taxon ID should be int')
    cursor.execute('SELECT `taxon_id` FROM taxon WHERE `ncbi_taxon_id`={};'
                   .format(starting_taxon))
    mysql_answer = cursor.fetchall()
    if not mysql_answer:
        raise ValueError('Invalid NCBI taxon id')
    mysql_id = mysql_answer[0][0]
    while True:
        cursor.execute('SELECT `parent_taxon_id` FROM taxon WHERE `taxon_id`={0};'
                       .format(mysql_id))
        mysql_answer = cursor.fetchall()
        mysql_id = mysql_answer[0][0]
        if mysql_id == 1:
            return
        else:
            cursor.execute('SELECT `ncbi_taxon_id` FROM taxon WHERE `taxon_id`={};'
                           .format(mysql_id))
            mysql_answer = cursor.fetchall()
            yield mysql_answer[0][0]

    
def get_taxa_list(taxon_id, cursor):
    """
    Return a list of taxa for a given taxon_id.
    Accepts a taxon id and returns a list of all taxa it belongs (as NCBI taxon
    ids), sorted from larger taxa to smaller ones. Does not include root
    :param taxon_id:
    :param cursor:
    :return:
    """
    #  Frankly, I'd like an excuse to cast it to list a third time in the same line
    return list(list(descend_taxon_tree(taxon_id, cursor)).__reversed__())


def is_taxon_member(taxon1, taxon2, cursor):
    """
    Return True if taxon 1 is a subtaxon of taxon2, False otherwise
    :param taxon_id:
    :param cursor:
    :return:
    """
    for taxon in descend_taxon_tree(taxon1, cursor):
        if taxon == taxon2:
            return True
    return False
    

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
    for supertaxon in descend_taxon_tree(taxon, cursor):
        if supertaxon in taxa_list:
            return supertaxon
    raise ValueError('Neither taxon is a supertaxon of a query')
