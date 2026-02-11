$infile =  $ARGV[0]; #snp50_file

$outfile = $ARGV[1];


open INFILE,  $infile;

open OUTFILE, ">$outfile";


$/ = "\n";

while  (<INFILE>)
{
my $line = $_;
chomp $line;
@tmp = split/\t/, $line;
#$ctg{$line} = $tmp[0];
$coord{$line} = $tmp[1];
push @contigs, $tmp[0];
push @lines, $line;
}

for $x (0..$#contigs)
{

if ($contigs[$x] eq $contigs[$x+1])
{
push @temp, "$coord{$lines[$x]}\t$lines[$x]";
}

else
{
push @temp, "$coord{$lines[$x]}\t$lines[$x]";
@temp = sort{$a<=>$b}@temp;

foreach $tmp (@temp)
{
print OUTFILE "$tmp\n";
}
@temp=();




}

}



