# Phylome

A collection of tools for phylomic analysis.

### Scripts

#### ~~clusterize.py~~ **OBSOLETE**
Takes a multiFASTA file(s) and does a simple distance-based clustering
using FFP distance metric.

#### cluster_diatoms.py
Clusters potentially homologous sequences based on BLAST output. The
algorithm is approximately similar to that of COG

#### filter_fasta.py
Takes a multiFASTA file(s) and removes sequences that are too short or
contain too many `X` characters.

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
