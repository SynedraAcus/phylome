A collection of tools for phylomic analysis. Depends on Biopython,
pymysql and numpy to work. As most of these scripts were written for a
single analysis, there may or may not be hardcoded values, questionable
interfaces and other such.

### Scripts

#### fasta_from_list.py
Takes FASTAs, the cluster file (every line should be a tab-separated
list of IDs) and BLAST file against and external DB and produces
separate FASTAs for every cluster and all its members' hits against that
DB

#### filter_fasta.py
Takes a multiFASTA file(s) and removes sequences that are too short or
contain too many `X` characters.

#### find_multiplicates.py
Takes a BLAST TSV file and (optionally) a multiFASTA. Returns IDs and
(if the FASTA is available) sequences of the multiplicated genes using
`phylome.multiplicates.is_duplicate`

#### parse_mapping.py
Processes IQtree likelihood mapping results. This script is meant to
integrate the results of multiple analyses into a single SVG image.
To make it feasible to process this image, separate data points are not
stored; instead, triangle is split into regions and quartet counts are
integrated over those.

#### process_outers.py
Iterates over clusters' hit fastas. Assembles statistics for
the amount of red or green (including higher plants) and splits the file
into red, green, and rest, filtering out diatoms.

#### one-shots/
A bunch of single-use scripts, mostly for reading results of something
or generating command lists. These *are not meant* for reuse.

### phylome

A collection of modules useful for a phylomic analysis.

#### phylome.blast_parser
A lightweight BLAST tabular output parser, written in a functional style
with loads of iterators. It is quicker than Biopython's BLAST modules,
but doesn't support formats other than tabular.

#### phylome.multiplicates
Detect intragenic duplications using BLAST hits against a database of
non-duplicated hits (probably usable with nr or any other huge reference
DB).

#### phylome.taxonomy
A few functions to query NCBI taxonomy database. Assumes database to be
stored in a mySQL instance
