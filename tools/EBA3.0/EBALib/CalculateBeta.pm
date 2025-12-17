# $Id$ Calculate Beta
# Perl module for EBA EBALib::CalculateBeta;
# Author: Jitendra Narayan <jnarayan81@gmail.com>
# Copyright (c) 2015 by Jitendra. All rights reserved.
# You may distribute this module under the same terms as Perl itself

##-------------------------------------------------------------------------##
## POD documentation - main docs before the code
##-------------------------------------------------------------------------##

=head1 NAME

EBALib::CalculateBeta - DESCRIPTION of Object

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

package EBALib::CalculateBeta;
use strict;
use warnings;
#use Term::ANSIColor;

use Exporter;

our @EXPORT_OK = "betaCal";

## Currently only the breakpoint graph is generated.

use EBALib::Draw::drawBreakpointGraph;
use EBALib::Draw::drawBetaGraph;
use EBALib::Draw::drawChrBreakpointGraph;
# use BreaksFinder; ## there is some changes in format so I decided to copy breafinder subrutine here manualy.

sub betaCal {

my ($dir, $increase, $num)=@_;

#print "Calculating Beta scores.. .\t\n";

if (-e "beta.tmp") { EBALib::Messages::betaHai();  unlink ("beta.tmp"); EBALib::Messages::betaDel();}   ## I check the directory and if the file already exist delete them .. becuase i am appending in file
if (-e "betaScore") { EBALib::Messages::betaScoreHai();  unlink ("betaScore"); EBALib::Messages::betaScoreDel();}   ## I check the directory and if the file already exist delete them .. becuase i am appending in file

my @SpsArray;
open SPSFILE, "sps.txt" or die $!;
while (<SPSFILE>) { my $SpsLine=$_; chomp $SpsLine; @SpsArray=split /,/, lc($SpsLine);  my $SpsNumber = scalar (@SpsArray); } ## It read the spacies names from sps.txt file ... need to improve !!!
close SPSFILE or die EBALib::Messages::failCl("sps.txt");

my $parent = "./$dir";  #The main folder which contains all the resolution folders
my ($par_dir, $sub_dir);
my @allResolutions;
opendir($par_dir, $parent) or die EBALib::Messages::failOp($parent); 
while (my $sub_folders = readdir($par_dir)) {
    next if ($sub_folders =~ /^\.\.?$/);  # skip . and ..
    my $path = $parent . '/' . $sub_folders;
    next unless (-d $path);   # skip anything that isn't a directory
	# my @name = split(/\//, $path);
	#print "Working with: $path .. .\t\n";
	generateBetaGancho($path, $increase);
	my @res = split(/\//, $path);
	my $resolutionName=$res[-1];  ## Name of the resolutions
	if (isInteger($resolutionName)) { 
		push @allResolutions, $resolutionName;
	}	
	else { 
		EBALib::Messages::numericFold($resolutionName);
	}
}

my @allResolutions_sorted = sort { $a <=> $b } @allResolutions; 

if (scalar(@allResolutions) <=2 ) { EBALib::Messages::noBeta();}      
# print @allResolutions_sorted;

###--- Calculate the breakpoint statics graph. -------

EBALib::Draw::drawBreakpointGraph::breakpointGraph ("beta.tmp", \@SpsArray, \@allResolutions, $num);
#EBALib::Draw::drawChrBreakpointGraph::breakpointChrGraph("beta.tmp", \@SpsArray, \@allResolutions);

##----------------------------------------------------

my $minResolution=$allResolutions_sorted[0];	
my $maxResolution=$allResolutions_sorted[-1];

# print "$minResolution\t$maxResolution\n";

for ( my $res=0; $res<=$#allResolutions_sorted; $res++) {
	if ($allResolutions_sorted[$res] == $minResolution) {
		lowerResCal($allResolutions_sorted[$res], $allResolutions_sorted[$res+1], $allResolutions_sorted[$res+2],\@SpsArray);
		}
	elsif ($allResolutions_sorted[$res] == $maxResolution) {
		#upperResCal($allResolutions_sorted[$res], $allResolutions_sorted[$res-1], $allResolutions_sorted[$res-2]);
		}
	else {
		middleResCal($allResolutions_sorted[$res-1], $allResolutions_sorted[$res], $allResolutions_sorted[$res+1], \@SpsArray);
		}
}

undef @allResolutions;

##------ Draw a graph for beta scores ----------

#EBALib::Draw::drawBetaGraph::drawBetaGraph;

##----------------------------------------------

unlink ("beta.tmp"); ## delete at the end  
closedir($par_dir);
}

## subroutines ===
##------------------------------------------------------------------------
sub middleResCal {
my ($minRes, $midRes, $maxRes, $spsArray_ref)=@_;
my @spsArray=@$spsArray_ref;
my (%allMinResHash, %allMidResHash, %allMaxResHash);
# print "$minRes, $midRes, $maxRes\n";
	foreach my $speciesName(@spsArray) {    # print "$speciesName\n";
	open INFILE,  'beta.tmp' or die EBALib::Messages::failOp('beta.tmp');
	open OUTFILE, ">>", 'betaScore' or die EBALib::Messages::failOp('betaScore');
	while (<INFILE>) {
		chomp;    
		my $line= $_;  #print "$line\n";
		$line=EBALib::CommonSubs::trim($line);
		my @tmp = split /\t/, $line;
		s{^\s+|\s+$}{}g foreach @tmp;
		next if $tmp[3] ne $speciesName;
		next if $tmp[6] ne "Break";
		if ($tmp[0] == $minRes) {
			my $minCoordinates="$tmp[4]<->$tmp[5]";
			$allMinResHash{$minCoordinates} = $tmp[2];
	        }
		elsif ($tmp[0] == $midRes) {
			my $midCoordinates="$tmp[4]<->$tmp[5]";
			$allMidResHash{$midCoordinates} = $tmp[2];
	        }
		elsif ($tmp[0] == $maxRes) {
			my $maxCoordinates="$tmp[4]<->$tmp[5]";
			$allMaxResHash{$maxCoordinates} = $tmp[2];
	        }
	}                   
        close INFILE;

my $realBreakpoints=compareResolutionsForAnyTwoMid(\%allMinResHash, \%allMidResHash, \%allMaxResHash, $midRes);
#print "$realBreakpoints\n";
my $missedBreakpoints=compareResolutionsForMissedMid(\%allMinResHash, \%allMidResHash, \%allMaxResHash, $midRes);
#print "$missedBreakpoints\n";
my $total=$missedBreakpoints+$realBreakpoints;
if($total ==0) { EBALib::Messages::zero(); $total = 1;}
my $finalBetaScore=$missedBreakpoints/$total;
# print "$missedBreakpoints/($missedBreakpoints+$realBreakpoints)\n";
#print "$finalBetaScore\n";
print OUTFILE "$midRes:$speciesName\t$finalBetaScore\n";
# print "$_ $allMidResHash{$_}\n" for (keys %allMidResHash);  ## to print hash
undef %allMinResHash; undef %allMidResHash; undef %allMaxResHash;
close OUTFILE; 
	}

}


##------------------------------------------------------------------------
sub lowerResCal {
my ($minRes, $midRes, $maxRes, $spsArray_ref)=@_;
my @spsArray=@$spsArray_ref;
my (%allMinResHash, %allMidResHash, %allMaxResHash);

	foreach my $speciesName(@spsArray) {    # print "$speciesName\n";
	open INFILE,  'beta.tmp' or die EBALib::Messages::failOp('beta.tmp');
	open OUTFILE, ">>", 'betaScore' or die EBALib::Messages::failOp('betaScore');
	while (<INFILE>) {
		chomp;    
		my $line= EBALib::CommonSubs::trim($_);  #print "$line\n";
		my @tmp = split /\t/, $line;
		s{^\s+|\s+$}{}g foreach @tmp;
		next if $tmp[3] ne $speciesName;
		next if $tmp[6] ne "Break";
		if ($tmp[0] == $minRes) {
			my $minCoordinates="$tmp[4]<->$tmp[5]";
			$allMinResHash{$minCoordinates} = $tmp[2];
	        }
		elsif ($tmp[0] == $midRes) {
			my $midCoordinates="$tmp[4]<->$tmp[5]";
			$allMidResHash{$midCoordinates} = $tmp[2];
	        }
		elsif ($tmp[0] == $maxRes) {
			my $maxCoordinates="$tmp[4]<->$tmp[5]";
			$allMaxResHash{$maxCoordinates} = $tmp[2];
	        }
	}
        close INFILE;

my $realBreakpoints=compareResolutionsForAnyTwo(\%allMinResHash, \%allMidResHash, \%allMaxResHash, $minRes);
#print "$realBreakpoints\n";
my $missedBreakpoints=compareResolutionsForMissed(\%allMinResHash, \%allMidResHash, \%allMaxResHash, $minRes);
#print "$missedBreakpoints\n";
my $total=$missedBreakpoints+$realBreakpoints;
if($total ==0) { EBALib::Messages::zero(); $total = 1;}
my $finalBetaScore=$missedBreakpoints/$total;
#print "$finalBetaScore\n";
print OUTFILE "$minRes:$speciesName\t$finalBetaScore\n";
# print "$_ $allMidResHash{$_}\n" for (keys %allMidResHash);  ## to print hash
undef %allMinResHash; undef %allMidResHash; undef %allMaxResHash;
close OUTFILE; 
	}

}

##----------------------------------------------------------------------------

sub compareResolutionsForMissed {

my ($allMinResHash_ref, $allMidResHash_ref, $allMaxResHash_ref, $minRes)=@_;
my %allMinResHash= %$allMinResHash_ref;
my %allMidResHash= %$allMidResHash_ref;
my %allMaxResHash= %$allMaxResHash_ref;
my $countMissed=0;
	# LOOP THROUGH IT
	while ( my ($key, $value) = each(%allMidResHash)){
     		#print $key.", ".$value."\n";
		my $resultMiss=checkOverlaps($value, $key, \%allMaxResHash);
			if (($resultMiss) == 1) { 
				my $finalResultMiss=checkOverlaps($value, $key, \%allMinResHash); 
				if (($finalResultMiss) == 0) { $countMissed++;}
				} 
	}
	
if ($countMissed > 0) { return $countMissed;} else { return 0;}
#print "$countReal\n";
undef %allMinResHash; undef %allMidResHash; undef %allMaxResHash;
}
##----------------------------------------------------------------------------

sub compareResolutionsForMissedMid {

my ($allMinResHash_ref, $allMidResHash_ref, $allMaxResHash_ref, $minRes)=@_;
my %allMinResHash= %$allMinResHash_ref;
my %allMidResHash= %$allMidResHash_ref;
my %allMaxResHash= %$allMaxResHash_ref;
my $countMissed=0;
	# LOOP THROUGH IT
	while ( my ($key, $value) = each(%allMinResHash)){
     		#print $key.", ".$value."\n";
		my $resultMiss=checkOverlaps($value, $key, \%allMaxResHash);
			if (($resultMiss) == 1) { 
				my $finalResultMiss=checkOverlaps($value, $key, \%allMidResHash); 
				if (($finalResultMiss) == 0) { $countMissed++;}
				} 
	}
	
if ($countMissed > 0) { return $countMissed;} else { return 0;}
#print "$countReal\n";
undef %allMinResHash; undef %allMidResHash; undef %allMaxResHash;
}

##----------------------------------------------------------------------------

sub compareResolutionsForAnyTwo {

my ($allMinResHash_ref, $allMidResHash_ref, $allMaxResHash_ref, $minRes)=@_;
my %allMinResHash= %$allMinResHash_ref;
my %allMidResHash= %$allMidResHash_ref;
my %allMaxResHash= %$allMaxResHash_ref;
my $countReal=0;
	# LOOP THROUGH IT
	while ( my ($key, $value) = each(%allMinResHash)){
     		#print $key.", ".$value."\n";
		my $resultMid=checkOverlaps($value, $key, \%allMidResHash);
		my $resultMax=checkOverlaps($value, $key, \%allMaxResHash);
		if (($resultMid || $resultMax) == 1) { $countReal++;}  
		}
return $countReal;
#print "$countReal\n";
undef %allMinResHash; undef %allMidResHash; undef %allMaxResHash;
}
##----------------------------------------------------------------------------

sub compareResolutionsForAnyTwoMid {

my ($allMinResHash_ref, $allMidResHash_ref, $allMaxResHash_ref, $midRes)=@_;
my %allMinResHash= %$allMinResHash_ref;
my %allMidResHash= %$allMidResHash_ref;
my %allMaxResHash= %$allMaxResHash_ref;
my $countReal=0;
	# LOOP THROUGH IT
	while ( my ($key, $value) = each(%allMidResHash)){
     		#print $key.", ".$value."\n";
		my ($resultMin, $numberMin)=checkOverlaps($value, $key, \%allMinResHash);
		my ($resultMax, $numberMax)=checkOverlaps($value, $key, \%allMaxResHash);
		if (($resultMin || $resultMax) == 1) { 
			$countReal++;
		}
	}
return $countReal;
#print "$countReal\n";
undef %allMinResHash; undef %allMidResHash; undef %allMaxResHash;
}

##------------------------------------------------------------------------------

sub checkOverlaps {
   #my ($brk,$brk_array_ref)= @_;
   my ($chr, $brkCoordinate, $hash_ref)=@_;
   my %hash=%$hash_ref;
   my @val1=split(/\<\-\>/, $brkCoordinate);	
   my @all_overlaps;  my $num=0;

	while ( my ($key, $value) = each %hash ){
  	#print "key: $key, value: $hash{$key}\n";  
	next if $value ne $chr;
	my @val2=split(/\<\-\>/, $key);
	my $OverRes = EBALib::CommonSubs::checkCorOverlaps ($val1[0],$val1[1],$val2[0],$val2[1]);
	if ($OverRes) {  
   		#push @all_overlaps, "$val2[0]\t$val2[1]";     
		$num++;
		}
	 }
if ($num>0) { return 1;} else { return 0;}
undef %hash; 
}

## ---------------------------------------------------------------------------- 
#subroutines here

sub isInteger { defined $_[0] && $_[0] =~ /^[+-]?\d+$/; }  ## It may accept - or + sign :(

## ----------------------------------------------------------------------------

sub generateBetaGancho {
my ($dir, $increase)=@_;
if(!$increase) { $increase=0;} ## The default values for increament is 0;
# print "Calculating Beta Scores .....\n";

		 opendir(DIR, $dir) or die $!;
   		 while (my $file = readdir(DIR)) {
        		# We only want files
        		next unless (-f "$dir/$file");
        		# Use a regular expression to find files ending in .txt
        		next unless ($file =~ m/\.txt$/); 
			my @new_file=split(/_/, $file);   ## file should be seperated by "_" !!!
				if ($new_file[0] ne "") { 
			        	# print "$file\n";
					my $in="$dir/$file";
					my $out="beta.tmp";
					#my $increase=0; ## It it the values accepted when user enter any option to increase
					findBreaks($dir,$in,$out,$increase);
				}
		  }
closedir(DIR);
}

##-------------------------------------------------------------------------------

sub findBreaks {

#!/usr/bin/perl
use strict;
use warnings;

my ($path, $InFile, $OutFile, $increase)=@_;

my (@array, @information, %scaffData);
## it sort the file with 2,3,4 columns and generate the breakpoints. !!!
## perl brk_finder.pl EBA_input_files/$file $new_file[0]"."_brk.eba      

my @res = split(/\//, $path);
my $resolutionName=$res[-1];  ## Name of the resolutions

open INFILE,  $InFile or die EBALib::Messages::failOp($InFile);
open OUTFILE, ">>" , $OutFile or die EBALib::Messages::failOp($OutFile); ##appending all the files in one tmp file

$|++;
$/ = "\n";
#my @array = <INFILE>; # Read all file to an array
while (<INFILE>) {
	my $line=lc($_); chomp $line;
	if ($line =~ /^\s*#/) { next; }
	$line=EBALib::CommonSubs::trim($line);
	my @tmpLine=split /\t/, $line;
	push (@information, $tmpLine[9]);  ## the scaffolds information or genome information
	push (@array, $line);
	$scaffData{"$tmpLine[5]:$tmpLine[6]"}=$tmpLine[4]; #Stored the name and coordinates information of scaffolds
}
close INFILE or die EBALib::Messages::failCl($InFile);

my @info=EBALib::CommonSubs::uniq(@information); 
if(scalar (@info) > 1) { EBALib::Messages::multiHit($InFile);}

my $status=lc($info[0]);
if ((scalar(@info) < 1) or (!@info)){ EBALib::Messages::genomeMSG($InFile); }

my @sorted_array = sort { (split "\t", $a)[1] cmp (split "\t", $b)[1] || (split "\t", $a)[2] <=> (split "\t", $b)[2] && (split "\t", $a)[3] <=> (split "\t", $b)[3] } @array;
push @sorted_array, "==="; ## End line to terninate;

if($#sorted_array < 2) { EBALib::Messages::fileFormat(); } ## If empty

# print map { "$_ $scaffData{$_}\n" } keys %scaffData;

	if ($status eq 'chromosomes') {
		createBreaks(\@sorted_array, $increase, $InFile, $resolutionName);
	}
	elsif ($status eq 'scaffolds') {
		createBreaksScaffolds(\@sorted_array, $increase, $InFile, \%scaffData, $resolutionName);
	}

undef @sorted_array; undef @array; undef @information; undef @info; undef %scaffData;
close OUTFILE or die EBALib::Messages::failCl($OutFile);

}  # subrutine ends here 


sub createBreaks {

my ($sorted_array_ref, $increase, $InFile, $resolutionName)= @_;
my @sorted_array=@$sorted_array_ref;

my $count; my $chr_num;
my (@ref_org, @chr_ref, @brk_cor1, @brk_cor2, @tar_org);


foreach my $line (@sorted_array) {
	chomp $line;
 	$line=lc($line);
	if ($line =~ /^\s*#/) { next; }
	$count++;
	my @tmp = split/\t/, lc($line);
	for (@tmp)  { s/^\s+//; s/\s+$//; } #replace one or more spaces 
	if ($count == 1) { $chr_num = $tmp[1];}

	##--------------------- incerase or decrease the breakpoints size start -----------
 
	if ($line =~ /^\=/) { $tmp[1]=0; $tmp[2]=0; $tmp[3]=0; }         ## change here in increments .. ????
	my $tmp3cor1=$tmp[3]-$increase;      
	my $tmp2cor2=$tmp[2]+$increase;

	##-------------------- incerase or decrease the breakpoints size end ----------- 

		if ($tmp[1] eq $chr_num) {
			push @ref_org, $tmp[0];
			push @chr_ref, $tmp[1];
			push @brk_cor1, $tmp3cor1;
			push @brk_cor2, $tmp2cor2;

			push @tar_org, $tmp[8];
			}

		else    {
			for my $x(0.. $#ref_org) {
				if ($x != $#ref_org) {
					print OUTFILE "$resolutionName\t$ref_org[$x]\t$chr_ref[$x]\t$tar_org[$x]\t$brk_cor1[$x]\t$brk_cor2[$x+1]\tBreak\n";
					# print  "$ref_org[$x]\t$chr_ref[$x]\t$tar_org[$x]\t$brk_cor1[$x]\t$brk_cor2[$x+1]\n";
						if ($brk_cor1[$x] > $brk_cor2[$x+1]) { 
							EBALib::Messages::confused($InFile,$chr_ref[$x],$brk_cor1[$x],$brk_cor2[$x+1]); 
							exit(1); ## Terminate the program
						}
					}
				}
			$chr_num=$tmp[1];

			undef @ref_org; undef @chr_ref; undef @brk_cor1; undef @brk_cor2; undef @tar_org;

			push @ref_org, $tmp[0];
			push @chr_ref, $tmp[1];
			push @brk_cor1, $tmp3cor1;
			push @brk_cor2, $tmp2cor2;

			push @tar_org, $tmp[8];

			}
}
undef @ref_org; undef @chr_ref; undef @brk_cor1; undef @brk_cor2; undef @tar_org;
} ##create breaks ends here 


sub createBreaksScaffolds {

my ($sorted_array_ref, $increase, $InFile, $scaffData_ref, $resolutionName)= @_;
my @sorted_array=@$sorted_array_ref;
my %scaffData=%$scaffData_ref;

#print map { "$_ $scaffData{$_}\n" } keys %scaffData;

my $count; my $chr_num;
my (@ref_org, @chr_ref, @brk_cor1, @brk_cor2, @tar_org, @scaffId, @scaffStart, @scaffEnd, @sign);

foreach my $line (@sorted_array) {
	chomp $line;
 	$line=lc($line);
	if ($line =~ /^\s*#/) { next; }
	$count++;
	my @tmp = split/\t/, lc($line);
	for (@tmp)  { s/^\s+//; s/\s+$//; } #replace one or more spaces 
	if ($count == 1) { $chr_num = $tmp[1];}

	##--------------------- incerase or decrease the breakpoints size start -----------
 
	if ($line =~ /^\=/) { $tmp[1]=0; $tmp[2]=0; $tmp[3]=0; }         ## change here in increments
	my $tmp3cor1=$tmp[3]-$increase;      
	my $tmp2cor2=$tmp[2]+$increase;

	##-------------------- incerase or decrease the breakpoints size end ----------- 

		if ($tmp[1] eq $chr_num) {
			push @ref_org, $tmp[0];
			push @chr_ref, $tmp[1];
			push @brk_cor1, $tmp3cor1;
			push @brk_cor2, $tmp2cor2;

			push @tar_org, $tmp[8];
			push @scaffId, $tmp[4];
			push @scaffStart, $tmp[5];
			push @scaffEnd, $tmp[6];
			push @sign, $tmp[7];
			}

		else    {
			for my $x(0.. $#ref_org) {
				if ($x != $#ref_org) {
					my $Id=$scaffId[$x];
					my @AllValuesFirst; my $decision; my $pahalaCor=0; my $dusaraCor=0; my $flagFirst;
					my @keysFirst = grep { $scaffData{$_} eq $Id } keys %scaffData; ## grep key if values match
						foreach (@keysFirst) { my @KvaluesFirst=split (/\:/, $_); foreach (@KvaluesFirst) { push @AllValuesFirst, $_;} }
						@AllValuesFirst=sort {$a <=> $b} (@AllValuesFirst);
						#foreach (@AllValuesFirst) { print "$_\t"; } print "\n"; exit;
						my @scaffBeginEndFirst = ($AllValuesFirst[0], $AllValuesFirst[-1]);
						if ($sign[$x] eq "-") { $pahalaCor=$scaffStart[$x];} else { $pahalaCor=$scaffEnd[$x];}

							my $resFirst=EBALib::CommonSubs::isInList($pahalaCor, @scaffBeginEndFirst); 
							if ((scalar(@AllValuesFirst) > 2) and ($resFirst == 0)) {
								$flagFirst="YES";
							}
							else {
								$flagFirst="NO";
							}
					my @AllValuesSecond; my $nextId=$x+1; my $flagSecond; my $newId =$scaffId[$nextId];
					my @keysSecond = grep { $scaffData{$_} eq $newId } keys %scaffData; ## grep key if values match
						foreach (@keysSecond) { my @KvaluesSecond=split (/\:/, $_); foreach (@KvaluesSecond) { push @AllValuesSecond, $_;} }
						@AllValuesSecond=sort {$a <=> $b} (@AllValuesSecond);
						
						my @scaffBeginEndSecond = ($AllValuesSecond[0], $AllValuesSecond[-1]);
						if ($sign[$nextId] eq "-") { $dusaraCor=$scaffEnd[$nextId];} else { $dusaraCor=$scaffStart[$nextId];}
						my $resSecond=EBALib::CommonSubs::isInList($dusaraCor, @scaffBeginEndSecond);

							if ((scalar(@AllValuesSecond) > 2) and ($resSecond == 0)) {
								$flagSecond="YES";
							}
							else {
								$flagSecond="NO";
							}
						#if ($sign[$nextId] eq "-") { $dusaraCor=$scaffEnd[$nextId];} else { $dusaraCor=$scaffStart[$nextId];}
						
							if (($flagFirst eq "YES") or ($flagSecond eq "YES")) { 
								$decision = "Break";
							} 
							else { 
								$decision = "PseudoBreak";
							}
							undef @AllValuesFirst;
							undef @AllValuesSecond;
							print OUTFILE "$resolutionName\t$ref_org[$x]\t$chr_ref[$x]\t$tar_org[$x]\t$brk_cor1[$x]\t$brk_cor2[$x+1]\t$decision\n";
								if ($brk_cor1[$x] > $brk_cor2[$x+1]) { 
									EBALib::Messages::confused($InFile,$chr_ref[$x],$brk_cor1[$x],$brk_cor2[$x+1]);
							exit(1); ## Terminate the program
						
						}
					
					}
				}
			$chr_num=$tmp[1];

			undef @ref_org; undef @chr_ref; undef @brk_cor1; undef @brk_cor2; undef @tar_org;undef @scaffId; undef @scaffStart; undef @scaffEnd; undef @sign;

			push @ref_org, $tmp[0];
			push @chr_ref, $tmp[1];
			push @brk_cor1, $tmp3cor1;
			push @brk_cor2, $tmp2cor2;

			push @tar_org, $tmp[8];
			push @scaffId, $tmp[4];
			push @scaffStart, $tmp[5];
			push @scaffEnd, $tmp[6];
			push @sign, $tmp[7];
			}
}
undef @ref_org; undef @chr_ref; undef @brk_cor1; undef @brk_cor2; undef @tar_org; undef @scaffId; undef @scaffStart; undef @scaffEnd; undef @sign;
} ##create breaks ends here 

1;
