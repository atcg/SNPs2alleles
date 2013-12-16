SNPs2alleles
============

This program takes diploid biallelic SNP genotype data and converts it into population-level allele count data for use in the program TreeMix (https://code.google.com/p/treemix/).

The input file should contain a header row with column names, and each row below that corresponds to one individual. The first column is the sample name, and the second column is the population name. All remaining rows are comma-separated genotypes for each locus. See the included data file (SNPdata.txt) for an example of how input data must be arranged.

Usage:
perl SNPs2alleles.pl --in <file> --pops <int> --loci <int> --out <file>

 Options:
   -in=s            Input filename
   -pops=i          Total number of populations
   -loci=i          Total number of loci
   -out=s           Output filename
   -help|man        Prints out documentation
