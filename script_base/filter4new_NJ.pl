#!/usr/bin/perl
use strict;
use warnings;

# Check arguments
if (@ARGV != 3) {
    die "Usage: perl script.pl <input_file> <output_file> <res_size>\n";
}

my ($infile, $outfile, $res_size) = @ARGV;

# Open files
open my $in,  '<', $infile  or die "Cannot open '$infile': $!\n";
open my $out, '>', $outfile or die "Cannot open '$outfile': $!\n";

# Arrays to store data
my (@lines, @coord, @status, @chr);
my %idd;   # coordinates to exclude

# Read input file
while (my $line = <$in>) {
    chomp $line;
    my @tmp = split /\t/, $line;

    push @lines,  $line;
    push @coord,  $tmp[2];  # coordinate
    push @status, $tmp[9];  # "BOUND" or "inter"
    push @chr,    $tmp[1];  # chromosome
}

close $in;

# Step 1: process "BOUND" blocks to mark coordinates to exclude
my @block_coords;
my $in_block = 0;

for my $x (0 .. $#coord) {
    if ($status[$x] eq "BOUND" && !$in_block) {
        @block_coords = ($coord[$x]);
        $in_block = 1;
    }
    elsif ($in_block && $status[$x] eq "inter") {
        push @block_coords, $coord[$x];
    }
    elsif ($in_block && $status[$x] eq "BOUND") {
        push @block_coords, $coord[$x];
        $in_block = 0;

        # sort and compute block size
        @block_coords = sort { $a <=> $b } @block_coords;
        my $min  = $block_coords[0];
        my $max  = $block_coords[-1];
        my $size = $max - $min;

        # mark coordinates if block <= 300000
        if ($size <= $res_size) {
            for my $c (@block_coords) {
                $idd{"$c,$chr[$x]"} = 1;
            }
        }

        @block_coords = ();
    }
}

# Step 2: check isolated BOUND triplets and mark
for my $x (1 .. $#coord - 1) {  # avoid edges
    if ($status[$x-1] eq "BOUND" && $status[$x] eq "BOUND" && $status[$x+1] eq "BOUND"
        && $chr[$x-1] eq $chr[$x] && $chr[$x] eq $chr[$x+1]) {

        $idd{"$coord[$x],$chr[$x]"} = 1;
    }
}

# Step 3: print remaining coordinates with index
my $i = 1;
for my $x (0 .. $#lines) {
    unless ($idd{"$coord[$x],$chr[$x]"}) {
        print $out "$lines[$x]\t$i\n";
        $i++;
    }
}

close $out;

print "Processing complete. Output written to '$outfile'.\n";
