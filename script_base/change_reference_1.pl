$infile =  $ARGV[0]; #snp50_file

$outfile = $ARGV[1];


open INFILE,  $infile;

open OUTFILE, ">$outfile";


$/ = "\n";

while  (<INFILE>)
{
my $line = $_;
$line =~s/\r\n/\n/;
@tmp = split /\t/, $line;
$tmp[4] =~s/Chr//;
$tmp[4] =~s/chr//;
print OUTFILE "$tmp[4]\t$tmp[5]\t$tmp[6]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$tmp[7]\t$tmp[8]\n";

}

