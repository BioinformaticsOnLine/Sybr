$infile =  $ARGV[0]; #snp50_file

$outfile = $ARGV[1];


open INFILE,  $infile;

open OUTFILE, ">$outfile";


$/ = "\n";

while  (<INFILE>)
{
my $line = $_;
chomp $line;
push @tmp, $line;
}

@tmp = sort @tmp;


foreach $tmp (@tmp)
{
print OUTFILE "$tmp\n"; 
}

