$infile =   $ARGV[0]; #hetero hap file

$outfile = $ARGV[1];
open INFILE,  $infile;
open OUTFILE,  ">$outfile";

$/ = "\n";

while  (<INFILE>)
{
$line = $_;
chomp $line;
@tmp = split /\t/, $line;
push @lines, $line;
push @coord, $tmp[5];
push @chr, $tmp[4];
push @chr1, $tmp[1];
#print $tmp[4];
}


for $x (0..$#coord)
{
if ((abs($coord[$x] - $coord[$x-1]) > 100000) or abs($coord[$x] - $coord[$x+1]) > 100000 or  ($chr[$x] ne $chr[$x-1]) or ($chr[$x] ne $chr[$x+1]) or ($chr1[$x] ne $chr1[$x-1]) or ($chr1[$x] ne $chr1[$x+1]))
{
print OUTFILE "$lines[$x]\tBOUND\n";
}
else
{
print OUTFILE "$lines[$x]\tinter\n";
}
}