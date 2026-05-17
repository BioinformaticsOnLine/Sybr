#!/usr/bin/perl
use strict;
use warnings;

# Check arguments
if (@ARGV != 2) {
    die "Usage: perl script.pl <input_file> <output_file>\n";
}

my ($infile, $outfile) = @ARGV;

# Open files safely
open my $in, '<', $infile or die "Cannot open input file '$infile': $!\n";
open my $out, '>', $outfile or die "Cannot open output file '$outfile': $!\n";

# Arrays to store data
my @records;

# Read input file
while (my $line = <$in>) {
    chomp $line;
    my @fields = split /\t/, $line;

    # Validate column count
    if (@fields < 9) {
        warn "Skipping line with fewer than 9 columns: $line\n";
        next;
    }

    # Store line as a hash
    push @records, {
        line  => $line,
        chr   => $fields[1],   # column 2
        start => $fields[5],   # column 6
        end   => $fields[6],   # column 7
        hap   => $fields[4],   # column 5
        sign  => $fields[8],   # column 9
    };
}

# Process lines
for my $i (0 .. $#records - 1) {
    my $current = $records[$i];
    my $next    = $records[$i + 1];

    # Case 1: Same chromosome, non-overlapping, '+' strand, same haplotype
    if ($current->{chr} eq $next->{chr} 
        && $current->{end} < $next->{start} 
        && $current->{sign} eq '+' 
        && $next->{sign} eq '+' 
        && $current->{hap} eq $next->{hap}) {
        print $out $current->{line}, "\n";
    }
    # Case 2: Same chromosome, overlapping, '-' strand, same haplotype
    elsif ($current->{chr} eq $next->{chr} 
        && $current->{end} > $next->{start} 
        && $current->{sign} eq '-' 
        && $next->{sign} eq '-' 
        && $current->{hap} eq $next->{hap}) {
        print $out $current->{line}, "\n";
    }
    # Case 3: Chromosome or strand changes
    elsif ($current->{chr} ne $next->{chr} 
        || $current->{sign} ne $next->{sign}) {
        print $out $current->{line}, "\n";
    }
}

# Always print the last line
print $out $records[-1]{line}, "\n";

# Close files
close $in;
close $out;

print "Processing complete. Output written to '$outfile'.\n";
