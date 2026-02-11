$infile =   $ARGV[0]; #hetero hap file

$outfile =   $ARGV[1]; #hetero hap file
open INFILE,  $infile;
open OUTFILE,  ">$outfile";

$/ = "\n";

while  (<INFILE>)
{
$line = $_;
chomp $line;
@tmp = split /\t/, $line;
push @lines, $line;
push @coord, $tmp[2]; #[5];
push @status, $tmp[9];
push @chr, $tmp[1];
}
$dd=0;
#print @status;
for $x (0..$#coord)
{
if ($status[$x] eq "BOUND" and $dd==0)
{
push @array, $coord[$x];
$dd=1;
#print "$coord[$x]\n";
}
elsif ($dd==1 and $status[$x] eq "inter")
{
push @array, $coord[$x];
#print "$coord[$x]\n";

} 
elsif ($dd==1 and $status[$x] eq "BOUND")
{
push @array, $coord[$x]; 
#print "$coord[$x]\n";
$dd=0; $s=0;
@array = sort{$a<=>$b}@array;
$min = $array[0];
$max = $array[$#array];
$size = $max-$min;
#print "$chr[$x] $min $max $size\n";
if ($size <=50000)
  {
  foreach $array (@array)
  {
   $idd{$array,$chr[$x]} =1;
   #print "$array\n";
  }

}


@array =();
#push @array, $coord[$x]; 
}
}

for $x (0..$#coord)
{
if ($status[$x-1] eq "BOUND" and $status[$x] eq "BOUND" and $status[$x+1] eq "BOUND" and $chr[$x-1] eq $chr[$x] and $chr[$x] eq $chr[$x+1])
{
$idd{$coord[$x],$chr[$x]} =1;
}
}


$i=1;
for $x (0..$#lines)
{
if ($idd{$coord[$x],$chr[$x]} != 1)
{
print OUTFILE "$lines[$x]\t$i\n";
$i++;
}
}