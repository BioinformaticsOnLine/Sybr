# $Id$ Breakpoint Scoring
# Perl module for EBA EBALib::BreaksScoring;
# Author: Jitendra Narayan <jnarayan81@gmail.com>
# Copyright (c) 2015 by Jitendra. All rights reserved.
# You may distribute this module under the same terms as Perl itself

##-------------------------------------------------------------------------##
## POD documentation - main docs before the code
##-------------------------------------------------------------------------##

=head1 NAME

EBALib::BreaksScoring  - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=cut

=head1 CONTACT

Jitendra <jnarayan81@gmail.com>

=head1 APPENDIX

The rest of the documentation details each of the object methods.

=cut

##-------------------------------------------------------------------------##
## Let the code begin...
##-------------------------------------------------------------------------##

package EBALib::BreaksScoring;
#use Term::ANSIColor;
use Exporter;

our @EXPORT_OK = "breakpointScoring1";

## This part of the scoring uses the old EBRs scoring approach

sub breakpointScoring1 {
my ($path, @fileNames)=@_;

my $dir="$path/EBA_OutFiles";
EBALib::Messages::callProb();

	opendir(DIR, $dir) or die $!;
   		while (my $file = readdir(DIR)) {
        	# We only want files
        	next unless (-f "$dir/$file");
        	# Use a regular expression to find files ending in .txt
        	next unless ($file =~ m/\.eba2$/); 
		my @new_file=split(/_/, $file);
		my @arg_score;
		if ($new_file[0] ne "") {
			if ( grep { $_ eq $new_file[0] } @fileNames ) { 
 				my @tmp_names=@fileNames;
				@tmp_names = grep { $_ ne $new_file[0] } @tmp_names;
				#map {delete $tmp_names[$_] if $_ eq $new_file[0] } @tmp_names; #it delete the desite values for array
				# @arg_score=("perl" ,"new_percent3.pl");
				my $mtfile="$path/EBA_OutFiles/"."$new_file[0]"."_table.eba2";
				push (@arg_score,$mtfile);
					foreach my $jitu (@tmp_names) {
						next if $jitu eq "";
						my $j="$path/EBA_OutFiles/"."$jitu"."_table.eba2";
						push (@arg_score,$j);
					}
        			push (@arg_score,$path);
				calculateScores(\@arg_score, $path);
				#system (@arg_score);
				# if ( $? == -1 ){ print "command failed: $!\n";}
				EBALib::Messages::build($new_file[0]);
			}
		}
	undef @arg_score;	
}
closedir(DIR);
}


sub calculateScores {
my ($InPutFiles_ref, $path)=@_;
my @InPutFiles=@$InPutFiles_ref;
my @resName=split /\//, $path;
my @t_name=split /\//, $InPutFiles[0]; @t_name=split /\./, $t_name[4];
my $file = "$InPutFiles[-1]/EBA_OutFiles/table_"."$t_name[0]".".eba3";
my $tmp = "$InPutFiles[-1]/EBA_OutFiles/table.tmp" . $$; # habitually I make tmp files unique in case of multiuser usage.

foreach my $argnum (0 .. $#InPutFiles-1) {

open  OUTFILE, ">$tmp" or die EBALib::Messages::failOp($tmp);
open INFILE, "$InPutFiles[$argnum]" or die EBALib::Messages::failOp($InPutFiles[$argnum]);

$|++;
$/="\n";

open SPSFILE, "sps.txt" or die $!;  ## It contain species names 
while (<SPSFILE>) {  $l= EBALib::CommonSubs::trim($_); chomp $l;  @t=split /,/, lc($l);  $t_len = scalar (@t); }
my $ts=join("\t", @t); 
close SPSFILE;

print OUTFILE  "Species\tChromosome\tBrk_Point\tBreakpoint_Decision\t$ts\t\tScore\n";
my ($spsName, $gaps, $breaks, $other, $telomere);
while (<INFILE>) {
	$line=$_;
	chomp $line;
	next if $line=~ m/^Species/; # Header of the file
	@tmp=split /\t/, $line;
	# @org=split /\_/,$tmp[0];
	$org[0] = $tmp[0];
	push (@chr, $tmp[1]);
	push (@brk_pt,$tmp[2]);
	push (@brk_decision, $tmp[3]);
	if ($tmp[3] eq "Breakpoints") { $breaks++;} elsif  ($tmp[3] eq "Gap") { $gaps++; } else { $other++; }
	my @corVal=split /\--/,$tmp[2]; if ($corVal[0] == 1) { $telomere++;}
	$spsName=$tmp[0];
	}
if ($argnum == 0) { #To check the data once for gaps and breaks
	my $output_file="gaps_brks.stats";
	my $pura=$breaks+$gaps+$other;
	open FH, ">>$output_file" or die EBALib::Messages::failOp($output_file);
	my $finalGap=$gaps-$telomere;
	if (-z "$output_file") { print FH "Resolution\tName\tTotal\tBreaksNum\tGapsNum\tTelomere\tfinalGap\tOther\n"; } #If file empty
	print FH "$resName[-1]\t$spsName\t$pura\t$breaks\t$gaps\t$telomere\t$finalGap\t$other\n";
	close FH or die EBALib::Messages::failCl("$output_file");
	}
close INFILE or die EBALib::Messages::failCl("InFile");

for (@brk_pt) { s/^\s+//; s/\s+$//;}  #replace one or more spaces at the end of i

if ($argnum <=0) { open INFILE2, $InPutFiles[0]; } else {open  INFILE2, "$file" or die EBALib::Messages::failOp("$file"); }

while  (<INFILE2>) {
	$flag=0; 
	$line1 = $_;
	chomp $line1;
	#print OUTFILE "$line1\n";
	@tmp1 = split /\t/, $line1;

	if ($line1=~ m/^Species/) {	
		$flag=1;
		foreach $xx(0..@tmp1) {
			if ("$tmp1[$xx]" eq "$org[0]") { $vv=$xx; }	
		}
	}
	@val_arr=split /\,/,$tmp1[$vv];

	for (@val_arr) { s/^\s+//; s/\s+$//; } #replace one or more spaces at the end of it
     	
	foreach $xyz(@val_arr)	{       
	#next if ($xyz == 0);
 	$f=0;
        	foreach $t1 (0..$#chr) {   
			if (($brk_pt[$t1] eq $xyz)  and  ($chr[$t1] eq $tmp1[1])) {
				# print "$brk_pt[$t1] eq $xyz)\n";
				#$tmp1[$vv]="";
				$tmp1[$vv]="$xyz"."="."$brk_decision[$t1]";
				push (@new_val, $tmp1[$vv]);
				#print OUTFILE "@tmp1\n";
				$f=1;
			}
		}
		if($f==0){ push (@new_val, $xyz); $f=0; }	
	}

## print OUTFILE "@new_val\n";  
@new_val=join (',',@new_val);	
 	if ($flag !=1){
		foreach (@tmp1) { 
			if ($cal == $vv){ print OUTFILE "@new_val\t";}
			else { print OUTFILE "$_\t";}
 			$cal++; 
			}
		print OUTFILE "\n";
		$cal=0;
		}
	undef @new_val;  undef @val_arr;
	}

close INFILE2 or die EBALib::Messages::failCl("InFile2");
close OUTFILE or die EBALib::Messages::failCl("OutFile");

rename ( $tmp,$file ); #love this function renaming   
undef @chr; undef @brk_pt; undef @brk_decision;    
#if ($argnum ==1) { exit;}
} 

} # Main subrutine ends here  
              
## Store the data

sub store {
my @row_data="";
my @val="";
@row_data=(@_);
foreach $q(@row_data){ if ($q =~ m/^0/){ @val=split /\,/, $q;} }
return @new_val;
}


1;


