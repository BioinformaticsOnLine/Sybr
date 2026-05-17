#!/usr/bin/perl
use strict;
use warnings;

# Check arguments
if (@ARGV != 3) {
    die "Usage: perl script.pl <input_file> <output_file> <res_size>\n";
}

my ($infile, $outfile, $res_size) = @ARGV;

# Open files safely
open my $in,  '<', $infile  or die "Cannot open '$infile': $!\n";
open my $out, '>', $outfile or die "Cannot open '$outfile': $!\n";

# Arrays to store lines, coordinates, and chromosomes
my (@lines, @coord, @chr);
my %idd;  # hash to track printed lines

# Read input file
while (my $line = <$in>) {
    chomp $line;
    my @tmp = split /\t/, $line;

    # Skip lines that do not have enough columns
    next unless @tmp >= 6;

    push @lines, $line;
    push @coord, $tmp[5];  # coordinate column (6th column)
    push @chr,   $tmp[1];  # chromosome column (2nd column)

    # Optional debug: print column 5 if needed
    # print "$tmp[4]\n";
}

# Process lines
for my $x (0 .. $#lines) {

    # Always print the first line
    if ($x == 0) {
        print $out "$lines[$x]\n";
        $idd{$lines[$x]} = 1;
        next;
    }

    # Print the line only if:
    # 1) Distance to previous line < $res_size
    # 2) Same chromosome
    # 3) Line not already printed
    if (abs($coord[$x] - $coord[$x-1]) <= $res_size
        && $chr[$x] eq $chr[$x-1]
        && !$idd{$lines[$x]}) {

        print $out "$lines[$x]\n";
        $idd{$lines[$x]} = 1;
    }

    # Lines that don't satisfy the above condition are skipped (no else)
}

close $in;
close $out;

print "Processing complete. Output written to '$outfile'.\n";
