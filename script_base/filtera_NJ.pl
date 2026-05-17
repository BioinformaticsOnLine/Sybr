#!/usr/bin/perl
use strict;
use warnings;

# Check for correct number of arguments
if (@ARGV != 2) {
    die "Usage: perl script.pl <input_file> <output_file>\n";
}

my ($infile, $outfile) = @ARGV;

# Open input and output files safely
open my $in, '<', $infile or die "Cannot open input file '$infile': $!\n";
open my $out, '>', $outfile or die "Cannot open output file '$outfile': $!\n";

# Arrays to store data
my @records;

# Read input file line by line
while (my $line = <$in>) {
    chomp $line;
    my @fields = split /\t/, $line;

    # Ensure there are enough columns in the input
    if (@fields < 5) {
        warn "Skipping line with fewer than 5 columns: $line\n";
        next;
    }

    # Store data in a hash for readability
    push @records, {
        line  => $line,
        chr   => $fields[1],  # column 2
        start => $fields[2],  # column 3
        end   => $fields[3],  # column 4
        hap   => $fields[4],  # column 5
    };
}

# Process records
for my $i (0 .. $#records - 1) {
    my $current = $records[$i];
    my $next    = $records[$i + 1];

    # Case 1: Same chromosome, non-overlapping, same haplotype
    if ($current->{chr} eq $next->{chr} 
        && $current->{end} < $next->{start} 
        && $current->{hap} eq $next->{hap}) {
        print $out $current->{line}, "\n";
    }
    # Case 2: Chromosome changes
    elsif ($current->{chr} ne $next->{chr}) {
        print $out $current->{line}, "\n";
    }
}

# Always print the last line
print $out $records[-1]{line}, "\n";

# Close filehandles
close $in;
close $out;

print "Processing complete. Output written to '$outfile'.\n";
