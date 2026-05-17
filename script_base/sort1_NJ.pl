#!/usr/bin/perl
use strict;
use warnings;

# Check arguments
if (@ARGV != 2) {
    die "Usage: perl script.pl <input_file> <output_file>\n";
}

my ($infile, $outfile) = @ARGV;

# Open files safely
open my $in,  '<', $infile  or die "Cannot open input file '$infile': $!\n";
open my $out, '>', $outfile or die "Cannot open output file '$outfile': $!\n";

# Data structures
my %by_contig;   # store lines grouped by contig
my @contig_order; # preserve contig order

while (my $line = <$in>) {
    chomp $line;
    my @cols = split /\t/, $line;

    my ($contig, $coord) = @cols[0,1];

    push @{ $by_contig{$contig} }, [$coord, $line];

    # keep track of the order of contigs
    push @contig_order, $contig unless exists $by_contig{$contig} && @{$by_contig{$contig}} > 1;
}

close $in;

# Process each contig
for my $contig (@contig_order) {
    my @lines_sorted = sort { $a->[0] <=> $b->[0] } @{ $by_contig{$contig} };

    for my $entry (@lines_sorted) {
        my ($coord, $line) = @$entry;
        print $out "$coord\t$line\n";
    }
}

close $out;

print "Processing complete. Output written to '$outfile'.\n";
