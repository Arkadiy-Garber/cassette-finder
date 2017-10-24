# gene-cassette-finder
Proteins that interact with each other in vivo are often transcribed together under the same operon. Consequently, these genes are often found in close proximity in genomes. This program determines whether a set of genes (provided in fasta format) are found encoded in close proximity in the genomes available on NCBI's nr database.

This program requires a Mac operating system, with Python, version 3 or later, installed.

Two arguments must be passed to this program:

    -One is a a file containing sequences of interest (i.e. sequences that you expect to be encoded closely together on the same chromosome) in fasta format.
  
    -The second argument that must be passed is the location of a database in which homologs to the genes will be searched for. This could be either the non-redundant (nr) NCBI database, or a folder containing genomes (in fasta format) in which you expect the gene cassettes to be encoded in. If your choice the latter, then you must the program with '-DB local'. Otherwise, the default option '-DB nr' will be called.
