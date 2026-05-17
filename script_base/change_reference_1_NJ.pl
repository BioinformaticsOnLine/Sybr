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

while (my $line = <$in>) {
    chomp $line;
    $line =~ s/\r$//;  # Remove carriage return if present

    my @tmp = split /\t/, $line;

    # Ensure we have enough columns
    next unless @tmp >= 9;

    # Remove 'Chr' or 'chr' from column 5 (index 4)
    $tmp[4] =~ s/^Chr//i;

    # Print selected columns to output
    print $out join("\t", $tmp[4], $tmp[5], $tmp[6], $tmp[1], $tmp[2], $tmp[3], $tmp[7], $tmp[8]), "\n";
}

close $in;
close $out;

print "Processing complete. Output written to '$outfile'.\n";
