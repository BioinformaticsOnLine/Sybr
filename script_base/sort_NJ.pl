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

# Read all lines into an array
my @lines = <$in>;
chomp @lines;          # Remove newline characters
@lines = sort @lines;  # Sort alphabetically

# Write sorted lines to output
foreach my $line (@lines) {
    print $out "$line\n";
}

close $in;
close $out;

print "Processing complete. Output written to '$outfile'.\n";
