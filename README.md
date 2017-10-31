# gene-cassette-finder
Proteins that interact with each other in vivo are often transcribed together under the same operon. Consequently, these genes are often found in close proximity in genomes. This program determines whether a set of genes (provided in FASTA format) are found encoded in close proximity in a given set of genomes.

# Dependencies
  -MacOS
  -Python3
  -BLAST

Usage:

You can provide this program with as many genes, in FASTA amino acid format, as you would like. These genes will be blasted across a database (for example, NCBI or RefSeq), and then clusetered according to acession numbers.
