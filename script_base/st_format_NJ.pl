#!/usr/bin/perl
use strict;
use warnings;

# Check arguments
if (@ARGV != 2) {
    die "Usage: perl script.pl <input_file> <output_file>\n";
}

my ($infile, $outfile) = @ARGV;

# Open files
open my $in,  '<', $infile  or die "Cannot open input file '$infile': $!\n";
open my $out, '>', $outfile or die "Cannot open output file '$outfile': $!\n";

my $i = 0;

while (my $line = <$in>) {
    chomp $line;
    my @tmp = split /\t/, $line;
    $i++;

    # Extract chr0 from $tmp[4]
    my $chr0;
    if ($tmp[4] =~ /(.+)_l/) {
        $chr0 = $1;
    }
    elsif ($tmp[4] =~ /scaffold_(.+)/) {
        $chr0 = $1;
    }
    else {
        $chr0 = $tmp[4];
    }

    # Extract chr from $tmp[1]
    my $chr;
    if ($tmp[1] =~ /scaffold_(.+)/) {
        $chr = $1;
    }
    elsif ($tmp[1] =~ /(.+)_/) {
        $chr = $1;
    }
    else {
        $chr = $tmp[1];
    }

    # Print formatted line
    print $out join("\t",
        $chr0,
        $tmp[5],
        $tmp[6],
        $i . $tmp[8],
        $i . $tmp[8],
        $i . $tmp[8],
        $chr,
        $tmp[2],
        $tmp[3]
    ), "\n";
}

close $in;
close $out;

print "Processing complete. Output written to '$outfile'.\n";
