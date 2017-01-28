This repository contains my solutions to three sequence alignment problems: pairwise, triple, and multiple sequence alignment. All three programs require the input sequences to be in the FASTA format and output the resulting alignment in the FASTA format. In all three algorithms, indel rates are converted to gap penalties by simply taking the log-10 of the specified indel rate. Also, I do not use affine gap penalties. The algorithms currently only use the BLOSUM62 matrix for (mis)match scores, but this may be expanded in the future.

* **[pairSeqAlign](pairSeqAlign.py): Perform pairwise sequence alignment**
    * This is simply an implementation of the Needleman-Wunsch algorithm (exact solution)
    * Usage: ``pairSeqAlign.py [-h] [-i INPUT] [-o OUTPUT] [-p INDEL]``
        * If no input file is specified (via `-i`), the program reads from standard input
        * If no output file is specified (via `-o`), the program writes to standard output
        * If no indel rate is specified, a default value of 0.1 is used (which corresponds to a gap penalty of -1)
    * Run ``./pairSeqAlign.py -h`` or ``./pairSeqAlign.py --help`` for usage help

* **[tripleSeqAlign](tripleSeqAlign.py): Perform triple sequence alignment**
    * This is simply a three-dimensional extension of the Needleman-Wunsch algorithm (exact solution)
    * Usage: ``tripleSeqAlign.py [-h] [-i INPUT] [-o OUTPUT] [-p INDEL]``
        * If no input file is specified (via `-i`), the program reads from standard input
        * If no output file is specified (via `-o`), the program writes to standard output
        * If no indel rate is specified, a default value of 0.1 is used (which corresponds to a gap penalty of -1)
    * Run ``./tripleSeqAlign.py -h`` or ``./tripleSeqAlign.py --help`` for usage help

* **[multiSeqAlign](multiSeqAlign.py): Perform multiple sequence alignment**
    * This is a progressive alignment heuristic (inexact solution)
    * First, a distance matrix is computed from the sequence data
        * If a tree is provided (via ``-t``), the distance matrix is computed from the tree
        * If a tree is not provided, the distance matrix is computed using the resulting scores of all n(n-1) pairwise sequence alignments
    * The initial alignment is constructed from the three closest sequences using [tripleSeqAlign](tripleSeqAlign.py)
    * One-by-one, additional sequences are added to the alignment in order of distance from the initial three sequences using a Needleman-Wunsch-like algorithm I created to align a single sequence with an existing multiple sequence alignment
    * Usage: ``multiSeqAlign.py [-h] [-i INPUT] [-o OUTPUT] [-p INDEL] [-t TREE]``
        * If no input file is specified (via `-i`), the program reads from standard input
        * If no output file is specified (via `-o`), the program writes to standard output
        * If no indel rate is specified, a default value of 0.1 is used (which corresponds to a gap penalty of -1)
    * Run ``./tripleSeqAlign.py -h`` or ``./tripleSeqAlign.py --help`` for usage help

REQUIREMENTS
===
Currently, both [pairSeqAlign](pairSeqAlign.py) and [tripleSeqAlign](tripleSeqAlign.py) are fully-functional without any additional dependencies. If no tree is specified in [multiSeqAlign](multiSeqAlign.py), it is also fully-functional without any additional dependencies, but if you wish to specify a tree file (which allows for a more efficient distance matrix computation), [DendroPy](http://www.dendropy.org/) must be installed.