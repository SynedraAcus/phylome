A collection of tools for phylomic analysis. Depends on Biopython and
numpy to work.

### Scripts

#### fasta_from_list.py
Takes FASTAs, the cluster file (every line should be a tab-separated
list of IDs) and BLAST file against and external DB and produces
separate FASTAs for every cluster and all its members' hits against that
DB

#### find_multiplicates.py
Takes a BLAST TSV file and (optionally) a multiFASTA. Returns IDs and
(if the FASTA is available) sequences of the multiplicated genes using
`phylome.multiplicates.is_duplicate`

#### filter_fasta.py
Takes a multiFASTA file(s) and removes sequences that are too short or
contain too many `X` characters.

#### process_outers.py
Iterates over clusters' hit fastas. Assembles statistics for
the amount of red or green (including higher plants) and splits the file
into red, green, and rest, filtering out diatoms.

### one-shots
A bunch of single-use scripts, mostly for reading one or another analysis
results. These are extremely quick-and-dirty and are unreliable at best
when used not for their original purpose.
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
