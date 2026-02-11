$infile =   $ARGV[0]; #hetero hap file

$outfile = $ARGV[1];

open INFILE,  $infile;

open OUTFILE, ">$outfile";

$/ = "\n";

while  (<INFILE>)
{
$line = $_;
chomp $line;
@tmp = split /\t/, $line;
push @lines, $line;
push @coord, $tmp[5];
push @chr, $tmp[1];
#print $tmp[4];
}


for $x (0..$#coord)
{
if ( (abs($coord[$x] - $coord[$x-1]) < 500000  and ($chr[$x] eq $chr[$x-1]) and $idd{$lines[$x]} !=1) or  (abs($coord[$x] - $coord[$x+1]) < 500000 and ($chr[$x] eq $chr[$x+1]) and $idd{$lines[$x]} !=1))

{

print OUTFILE "$lines[$x]\n";
$idd{$lines[$x]} =1;

}

}
