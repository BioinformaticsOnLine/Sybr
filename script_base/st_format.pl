$infile =   $ARGV[0]; #hetero hap file
open INFILE,  $infile;
$outfile =   $ARGV[1]; #hetero hap file
open OUTFILE,  ">$outfile";


$i=0;
$/ = "\n";

while  (<INFILE>)
{
$line = $_;
chomp $line;
@tmp = split /\t/, $line;
$i++;


if ($tmp[4] =~/\_/)
{ 
($chr0) = $tmp[4] =~/(.+)\_l/;
}
if ($tmp[4] =~/scaffold_.+/)
{ 
($chr0) = $tmp[4] =~/scaffold_(.+)/;
}

else
{
$chr0=$tmp[4];
}

if ($tmp[1] =~/scaffold_.+/)
{ 
($chr) = $tmp[1] =~/scaffold_(.+)/;
}

elsif ($tmp[1] =~/\_/)
{ 
($chr) = $tmp[1] =~/(.+)\_/;
}
else
{
$chr=$tmp[1];
}
print OUTFILE "$chr0\t$tmp[5]\t$tmp[6]\t$i$tmp[8]\t$i$tmp[8]\t$i$tmp[8]\t$chr\t$tmp[2]\t$tmp[3]\n";

}


