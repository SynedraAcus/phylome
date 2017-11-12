# Phylome

A collection of tools for phylomic analysis. Depends on Biopython and
numpy to work.

### Scripts

#### find_multiplicates.py
Takes a BLAST TSV file and (optionally) a multiFASTA. Returns IDs and
(if the FASTA is available) sequences of the multiplicated genes using
`phylome.multiplicates.is_duplicate`

#### filter_fasta.py
Takes a multiFASTA file(s) and removes sequences that are too short or
contain too many `X` characters.

### one-shots
A bunch of single-use scripts, mostly for reading one or another analysis
results. These are extremely quick-and-dirty and are unreliable at best
when used not for the original analysis.
### phylome

A collection of modules useful for a phylomic analysis.

#### phylome.blast_parser
A lightweight BLAST tabular output parser, written in a functional style
with loads of iterators. It is quicker than Biopython's BLAST modules,
but it doesn't support formats other than tabular.

#### phylome.multiplicates
Detect intragenic duplications using BLAST hits against a database of
non-duplicated hits (probably usable with nr or any other huge reference
DB).

#### phylome.taxonomy
A few functions to query NCBI taxonomy database. Assumes database to be
stored in a mySQL instance

#### mcl_clustering
A copy of [koteth's implementation of MCL](https://github.com/koteth/python_mcl)
Requires numpy to work.
