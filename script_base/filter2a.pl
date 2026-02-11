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
push @start, $tmp[5];
push @end,   $tmp[6];
push @chr,   $tmp[1];
push @chr1,  $tmp[4];
push @sign,  $tmp[8];
}


for $x (0..$#chr)
{
  if ($chr[$x] eq $chr[$x+1] and $end[$x] < $start[$x+1] and $sign[$x] eq "+" and $sign[$x+1] eq "+" and $chr1[$x] eq $chr1[$x+1])
    {
     print OUTFILE "$lines[$x]\n";
    }
  elsif ($chr[$x] eq $chr[$x+1] and $end[$x] > $start[$x+1] and $sign[$x] eq "-" and $sign[$x+1] eq "-" and $chr1[$x] eq $chr1[$x+1])
    {
     print OUTFILE "$lines[$x]\n";
    }

   elsif ($chr[$x] ne $chr[$x+1] or $sign[$x] ne $sign[$x+1])
    {  
     print OUTFILE "$lines[$x]\n";
    }
}
