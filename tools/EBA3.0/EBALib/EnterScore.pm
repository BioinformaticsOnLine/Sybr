# $Id$ Enter Score
# Perl module for EBA EBALib::EnterScore;
# Author: Jitendra Narayan <jnarayan81@gmail.com>
# Copyright (c) 2015 by Jitendra. All rights reserved.
# You may distribute this module under the same terms as Perl itself

##-------------------------------------------------------------------------##
## POD documentation - main docs before the code
##-------------------------------------------------------------------------##

=head1 NAME

EBALib::EnterScore  - DESCRIPTION of Object

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

package EBALib::EnterScore;
#use Term::ANSIColor;
use Exporter;

our @EXPORT_OK = "enterScore";

## Old EBRs calculation approach

sub enterScore {
my ($path, @fileNames)=@_;
my $dir="$path/EBA_OutFiles";
EBALib::Messages::classProb();

		 opendir(DIR, $dir) or die $!;
   		 while (my $file = readdir(DIR)) {
        		# We only want files
        		next unless (-f "$dir/$file");
        		# Use a regular expression to find files ending in .txt
        		next unless ($file =~ m/\.eba4$/); 
			my @new_file=split(/_/, $file);
			my @arg_score;
				if ($new_file[0] ne "") { 
					if ( grep { $_ eq $new_file[0] } @fileNames ) 
					{ 
					my @tmp_names=@fileNames;
					@tmp_names = grep { $_ ne $new_file[0] } @tmp_names;
					###map {delete $tmp_names[$_] if $_ eq $new_file[0] } @tmp_names; #it delete the desite values for array
					#@arg_score=("perl" ,"enter_score.pl");
					my $mtfile="$path/EBA_OutFiles/"."$new_file[0]"."_table_scored.eba4";
					push (@arg_score,$mtfile);
						foreach my $jitu (@tmp_names)
							{
							next if $jitu eq "";
							my $j="$path/EBA_OutFiles/"."$jitu"."_table_scored.eba4";
							push (@arg_score,$j);
							}
					push (@arg_score,$path);
					addScore (@arg_score);
        				#system (@arg_score);
					#if ( $? == -1 ){ print "command failed: $!\n";}
					}
				}
		  undef @arg_score;				
		  }

closedir(DIR);
}

sub addScore { 

my @InputSpeciesNames=@_;

#!/usr/bin/perl
#use strict;
#use warnings;

@t_name=split /\//, $InputSpeciesNames[0]; @t_name=split /\./, $t_name[4];
my $file = "$InputSpeciesNames[-1]/EBA_OutFiles/table_"."$t_name[0]"."_scored.eba5";

#my $file = "table_table_$ARGV[0].txt";
my $tmp = "$InputSpeciesNames[-1]/EBA_OutFiles/table.tmp" . $$; # habitually I make tmp files unique in case of multiuser usage.

foreach $argnum (0 .. $#InputSpeciesNames-1) {
 
 open  OUTFILE, ">$tmp";
 open INFILE, $InputSpeciesNames[$argnum];

$/ = "\n";

open SPSFILE, "sps.txt" or die $!;
while (<SPSFILE>) { $l=$_; chomp $l; @t=split /,/, lc($l); $t_len = scalar (@t); }
$ts=join("\t", @t); 
close SPSFILE or die EBALib::Messages::failCl("sps.txt");

print OUTFILE  "Species\tChromosome\tBrk_Point\tBreakpoint_Decision\t$ts\t\tScore\n";

my @org;

while (<INFILE>)
{
$line=$_;
chomp $line;
next if $line=~ m/^Species/;
@tmp=split /\t/, $line;
$org[0]=$tmp[0];    ### !!!!!!! changed here from spliting tmp[0] to store in @arg ...
push (@chr, $tmp[1]);
push (@brk_pt,$tmp[2]);
push (@brk_decision, $tmp[3]);
push (@score, $tmp[-1]); #### We need to change this values as the colunm always changes when we add species. !!!!!!!
}
close INFILE or die EBALib::Messages::failCl("file");

for (@brk_pt) { s/^\s+//; s/\s+$//;}

if ($argnum <=0) { open INFILE2, $InputSpeciesNames[0] or die EBALib::Messages::failOp($InputSpeciesNames[0]); } else {open  INFILE2, "$file" or die  EBALib::Messages::failOp($file); }

#open INFILE2, $ARGV[1];
while  (<INFILE2>) {
$flag=0;
$line1 = $_;
chomp $line1;
@tmp1 = split /\t/, $line1;

	if ($line1=~ m/^Species/) {	
		$flag=1;
		foreach $xx(0..@tmp1) {
			if ("$tmp1[$xx]" eq "$org[0]")	{
	 			$vv=$xx; 
				# print "$vv\n";
			}
		}
	}
        @val_arr1=split /\,/,$tmp1[$vv];

	for (@val_arr1) { s/^\s+//; s/\s+$//;} 

	#print OUTFILE "@val_arr1\n";
	foreach (@val_arr1){
		@val_arr2=split /\=/,$_;
		push (@val_arr,$val_arr2[0]);
	}

for (@val_arr) { s/^\s+//; s/\s+$//;}

#	print OUTFILE "@val_arr\n";
   		
	foreach $xyz(@val_arr){       
			#next if ($xyz == 0);
	  		$f=0;
                        foreach $t1 (0..$#chr) {
				if (("$brk_pt[$t1]" eq "$xyz" )  and  ("$chr[$t1]" eq "$tmp1[1]")) {
				#$tmp1[$vv]="";
				$tmp1[$vv]="$xyz"."="."$brk_decision[$t1]"."+"."$score[$t1]";
				push (@new_val, $tmp1[$vv]);
				# push (@new_val, $score[$t1]);
				#print OUTFILE "@tmp1\n";
				$f=1;
				}
			}
			if($f==0) { push (@new_val, $xyz); $f=0; }		
		}


  #print OUTFILE "@new_val\n";   	
@new_val=join (',',@new_val);	
 
if ($flag !=1) {
foreach (@tmp1) { 
	if ($cal == $vv) { print OUTFILE "@new_val\t"; }
	else { print OUTFILE "$_\t"; }
 	$cal++; 
	}
print OUTFILE "\n";
$cal=0;
}
 
undef @new_val;  undef @val_arr;     undef @val_arr1; undef @val_arr2;
}         
    
#print "writing----------------------------";
close INFILE2 or die EBALib::Messages::failCl("outfile1");
close OUTFILE or die EBALib::Messages::failCl("outfile2");

rename ( $tmp,$file ); #love this function renaming   
undef  @chr; undef @brk_pt; undef @brk_decision; undef @score;   
#if ($argnum ==1) { exit;}
} 

} # Subrutone ends herer  

 			

1;
	 
