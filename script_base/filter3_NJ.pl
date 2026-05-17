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

# Arrays to store data
my (@lines, @coord, @chr, @chr1);

# Read input
while (my $line = <$in>) {
    chomp $line;
    my @tmp = split /\t/, $line;

    # Skip lines that don’t have enough columns
    next unless @tmp >= 6;

    push @lines, $line;
    push @coord, $tmp[5];   # coordinate column
    push @chr,   $tmp[4];   # column 5
    push @chr1,  $tmp[1];   # column 2
}

# Process lines
for my $x (0 .. $#lines) {

    # Compute distances safely
    my $dist_prev = ($x > 0) ? abs($coord[$x] - $coord[$x-1]) : undef;
    my $dist_next = ($x < $#lines) ? abs($coord[$x] - $coord[$x+1]) : undef;

    # Determine if the line is a boundary
    if ((!defined $dist_prev || $dist_prev > $res_size)
        || (!defined $dist_next || $dist_next > $res_size)
        || ($chr[$x] ne $chr[$x-1] || $chr1[$x] ne $chr1[$x-1])
        || ($chr[$x] ne $chr[$x+1] || $chr1[$x] ne $chr1[$x+1])) {

        print $out "$lines[$x]\tBOUND\n";
    }
    else {
        print $out "$lines[$x]\tinter\n";
    }
}

close $in;
close $out;

print "Processing complete. Output written to '$outfile'.\n";
