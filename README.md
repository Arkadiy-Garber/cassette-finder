# gene-cassette-finder
In bacteria, proteins that interact with each other in vivo are often encoded under the same operon. Thus, bacterial genomes consist of groups of genetic nieghborhoods, or gene cassettes, that encode genes with related functions (e.g. energy metabolism, motility, communication). In this way, cells can more efficiently regulate multiple genes at once.

In addition to sequence homology between individual genes, the co-occurence of certain genes within particular gene cassettes may be another indicator of function within the cell. To allow for the detection of possible gene cassettes, this 'cassette\_finderV1.py' program searches (via BLAST) for genes of interest, and utilizes the sequential assession numbers of genes within each genome to cluster matches, and identify which homologs appear to be encoded closely to one another. You can also blast against a set of locally downloaded genomes, in which case, fake sequential assession numbers will be created for clustering.

To start, you will need a set of genes (or one gene, if you are interested in identifying multiple copies of the same gene within one genomic neightborhood) in FASTA format, and a BLAST search engine installed. The script also requires Python version 3, and utilizes OSX commands, so if you are running on anything besides a MAC, you will likely have some hurdles to overcome.

The database against which your which genes of interest will be blasted is also up to you. You can blast against a set of genomes that you have downloaded; in this case, please provide a folder that contains only the genomes (in FASTA format) in which you wish to search for gene cassettes. You can also search against the non-redundant NCBI database (nr) or RefSeq, but these databases would have to be locally downloaded.

# Dependencies
  -MacOS
  
  -Python3
  
  -BLAST

# Usage:


