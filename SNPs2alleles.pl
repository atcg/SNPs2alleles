#!/usr/bin/perl


# The MIT License (MIT)
# 
# Copyright (c) 2013 Evan Melstad (evanmelstad@ucla.edu)
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
# the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my $help = 0;
my $inFile;
my $outFile;
my $popNumber;
my $totalLoci;

GetOptions  ("in=s"      => \$inFile,
             "out=s"     => \$outFile,
             "pops=i"    => \$popNumber,
             "loci=i"    => \$totalLoci);

if (!$inFile or !$outFile or !$popNumber or !$totalLoci or $help) {
    pod2usage(-exitval => 0, -verbose => 2, -noperldoc => 1);
}


open(my $in, "<", $inFile);

# Create a hash of hashes with the following structure:
#   %indHash {
#       0 {
#           sample => HBS39823
#           pop => pop1
#           genotypes => [TT,CC,CC,CC,CT...]
#         }
#       1 {
#           sample => HBS121940
#           pop => pop1
#           genotypes => [TT,CT,CC,CC,TT...]
my %indHash;
my $counter = 0;
while (my $line = <$in>) {
    next if ($line =~ /^SAMPLE#\tPOP\t/);
    my @columns = split ("\t", $line);
    $indHash{$counter}{"sample"} = shift(@columns);
    $indHash{$counter}{"pop"} = shift(@columns);
    $indHash{$counter}{"genotypes"} = (\@columns); #An array where each value is a comma-separated list of alleles for each locus
    $counter++;

}

# We want to get:
#popHash
#   Pop1 {
#      L1 {
#          A => 4
#          T => 0
#          }
#      L2 {
#          C => 3
#          T => 1
#          }
#    }
#    Pop2 {
#      L1 {
#          A => 4
#          T => 2
#          }
#      L2 {
#          C => 3
#          T => 4
#          }

my %popHash;
my $counter2 = 0;
foreach my $ind (sort (keys %indHash)) { #iterate through every row of dataset
    my $indPop = $indHash{$ind}{"pop"}; # the population this individual belongs to
    
    # Add the population to %popHash with 0's for all allele counts for all loci unless it has already been added already
    unless (exists $popHash{$indPop}) {
        my $locCounter = 0;
        while ($locCounter < $totalLoci) {
            # $popHash{1}{0} is the allele counts for population 1 at locus 0 (column three on the data file). $popHash{1}{1} is the allele counts for population 1 at locus 1 (column four on the data file).
            $popHash{$indPop}{$locCounter} = {"A"=>0, "T"=>0, "C"=>0, "G"=>0};
            $locCounter++;
        }
    }
    my $locCounter2 = 0; # Restart the counter at locus 0 (FOXG1)
    foreach my $genoString (@{$indHash{$ind}{"genotypes"}}) { # @{$indHash{$ind}{"genotypes"}} is an array holding all the genotype values for $ind. $genoString is iterating through these elements
        chomp($genoString);
        my @alleles = split(",", $genoString); # Now alleles is an array with two values--one for each allele in an individual
        foreach my $allele (@alleles) {
            $popHash{$indPop}{$locCounter2}{$allele}++; # For each of the alleles the individual hash, incremenent the population allele counter for that locus
        }
        $locCounter2++;
    }

}
# print Dumper(\%popHash); # Turn this back on if things aren't looking right to see the structure of the final hash
print "Data processing complete. Creating output file.\n";


#Now create output file
open(my $outFH, ">", $outFile) or die "couldn't create $outFile: $!\n";
    
foreach my $populationz (sort { substr($a, 3) <=> substr($b, 3)} (keys %popHash)) { #{ substr($a, 3) <=> substr($b, 3)} is required to properly order populations labeled with "pop1, pop2, ... pop11, pop12" so it doesn't order them "pop1, pop10, pop11, ..., pop2, pop20..."
    print $outFH "$populationz\t"; # Print a header row of all the population names
}
print $outFH "\n";


my $locCounter3 = 0; #Each locus will get its own line now
while ($locCounter3 < $totalLoci) {
    foreach my $populationz (sort { substr($a, 3) <=> substr($b, 3) } (keys %popHash) ) { # The populations each get their own column
        my $successCount = 0; # Since there are only two alleles per locus, we need to keep track of how many non-zero counts there are to get the formatting right
        if ($popHash{$populationz}{$locCounter3}{"T"} != 0) {
            print $outFH $popHash{$populationz}{$locCounter3}{"T"} . ","; # This corresponds to the number of "T" alleles that were found for that locus in that population
            $successCount++;
        }
        if ($popHash{$populationz}{$locCounter3}{"C"} != 0) {
            print $outFH $popHash{$populationz}{$locCounter3}{"C"}; # This corresponds to the number of "C" alleles that were found for that locus in that population
            $successCount++;
            if ($successCount == 2) {
                print $outFH "\t";
            } else {
                print $outFH ",";
            }
        }    
        if ($popHash{$populationz}{$locCounter3}{"A"} != 0) {
            print $outFH $popHash{$populationz}{$locCounter3}{"A"}; # This corresponds to the number of "A" alleles that were found for that locus in that population
            $successCount++;
            if ($successCount == 2) {
                print $outFH "\t";
            } else {
                print $outFH ",";
            }        
        }    
        if ($popHash{$populationz}{$locCounter3}{"G"} != 0) {
            print $outFH $popHash{$populationz}{$locCounter3}{"G"}; # This corresponds to the number of "G" alleles that were found for that locus in that population
            $successCount++;
            if ($successCount == 2) {
                print $outFH "\t";
            } else {
                print $outFH ",";
            }        
        }
        if ($successCount == 0) {
            print $outFH "0,0\t"; # Required in case there is all missing data for a population at a locus
        } elsif ($successCount == 1) {
            print $outFH "0\t"; # Required in case a population is fixed for an allele at a locus
        }
    }
    print $outFH "\n";
    $locCounter3++;
}
print "\nFinished printing output file.\n";




#Documentation
__END__

=head1 NAME

SNPs2alleles.pl

=head1 SYNOPSIS 

perl SNPs2alleles.pl --in <file> --pops <int> --loci <int> --out <file>

 Options:
   -in=s            Input filename
   -pops=i          Total number of populations
   -loci=i          Total number of loci
   -out=s           Output filename
   -help|man        Prints out documentation


=head1 DESCRIPTION

This program takes diploid biallelic SNP genotype data and converts it into population-level allele count data for use in the program TreeMix (https://code.google.com/p/treemix/).

The input file should contain a header row with column names, and each row below that corresponds to one individual. The first column is the sample name, and the second column is the population name. All remaining rows are comma-separated genotypes for each locus. See the included data file (SNPdata.txt) for an example of how input data must be arranged.

=cut






