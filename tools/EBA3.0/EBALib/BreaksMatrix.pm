# $Id$ Breakpoint Matrix
# Perl module for EBA EBALib::BreaksMatrix;
# Author: Jitendra Narayan <jnarayan81@gmail.com>
# Copyright (c) 2015 by Jitendra. All rights reserved.
# You may distribute this module under the same terms as Perl itself

##-------------------------------------------------------------------------##
## POD documentation - main docs before the code
##-------------------------------------------------------------------------##

=head1 NAME

EBALib::BreaksMatrix  - DESCRIPTION of Object

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

package EBALib::BreaksMatrix;
use strict;
use warnings;
#use Term::ANSIColor;

use Exporter;

our @EXPORT_OK = "breakpointsTable";

sub breakpointsTable {
my $path=shift;
my $dir="$path/EBA_OutFiles";
EBALib::Messages::tableMSG();
		opendir(DIR, $dir) or die $!;
   		while (my $file = readdir(DIR)) {
        		# We only want files
        		next unless (-f "$dir/$file");
        		# Use a regular expression to find files ending in .txt
        		next unless ($file =~ m/\.eba1$/); 
			my @new_file=split(/_/, $file);
				if ($new_file[0] ne "") { 
					my $InFileName= "$path/EBA_OutFiles/$new_file[0]_brk.eba";
					my $AllFileName= "$path/EBA_OutFiles/$file";
					my $OutFileName= "$path/EBA_OutFiles/$new_file[0]_table.eba2";
					drawMatrix ($InFileName, $AllFileName, $OutFileName); 

				}				
		 }
		closedir(DIR);
}

sub drawMatrix {

my ($InFile, $InFile2, $OutFile)=@_;

open INFILE,  $InFile or die EBALib::Messages::failOp($InFile);
open OUTFILE, ">" , $OutFile or die EBALib::Messages::failOp($OutFile);

open SPSFILE, "sps.txt" or die $!;
my @SpsArray;
while (<SPSFILE>) { my $SpsLine=$_; chomp $SpsLine; @SpsArray=split /,/, lc($SpsLine);  my $SpsNumber = scalar (@SpsArray); } ## It read the spacies names from sps.txt file ... need to improve !!!
my $SpsArrayTabed=join("\t", @SpsArray); 
close SPSFILE or die EBALib::Messages::failCl("sps.txt");

print OUTFILE "Species\tChromosome\tBreakpoint_Coordinates\tBreakpoint_Decision\t$SpsArrayTabed\t\tScore\n";

$|++;
$/ = "\n";

my @matrix;
while  (<INFILE>) {   ## Read *.eba1 file for one species at a time
my $line = $_; chomp $line;
$line=EBALib::CommonSubs::trim($line);
my @tmp = split /\t/, $line;
my $lineNumber=$.;
my @allLines;
	open INFILE2,  $InFile2 or die EBALib::Messages::failOp($InFile2);                            
	while  (<INFILE2>) {   ## Read *.eba1 file for one species at a time
		my $line1 = $_; chomp $line1;
		$line1=EBALib::CommonSubs::trim($line1);
		my @tmp1 = split /\t/, $line1;
			if (($tmp[1] eq $tmp1[0]) and ($tmp[3] == $tmp1[3])) { push @allLines, $line1;  } 

		} ## end INFILE2 loop here
	close INFILE2 or die EBALib::Messages::failCl("file");

my @allUniqueLines=EBALib::CommonSubs::sortUniqueHash(@allLines);
my @matrix=enterValue(\@allUniqueLines, \@SpsArray, $lineNumber, $tmp[5]);

foreach my $value (@matrix) { print OUTFILE "$value\t";} print OUTFILE "\n"; 

undef @matrix; undef @allLines;	
}
close INFILE or die EBALib::Messages::failCl($InFile);
close OUTFILE or die EBALib::Messages::failCl($OutFile);

} ## the Main subrutines ends here


## all subrutines
##------------------------------------------------------------------------------------------
sub enterValue {
  my ($allUniq_ref, $SpsArray_ref, $lineNumber, $brkInfo)=@_;
  my  @allUniq=@$allUniq_ref; my @matrix;  my %doneHash; my %doneHashPseudo;
  my  @SpsArray=@$SpsArray_ref;  my $multipleHit="NO"; my @allCoordinates;
	for (my $i=0; $i<=$#SpsArray+4; $i++) { $matrix[$i]=0; } # Inserted zeros
	foreach my $line (@allUniq) {
		$line=EBALib::CommonSubs::trim($line);
		my @tmp = split /\t/, $line;

#1	pygoscelis_adeliae	picoides_pubescens	588086	671505	615977	627480	PseudoBreak	PseudoBreak
#1	pygoscelis_adeliae	egretta_garzetta	588086	671505	423266	658440	PseudoBreak	PseudoBreak
#1	pygoscelis_adeliae	manacus_vitellinus	588086	671505	574999	658406	PseudoBreak	Break

		# my @tarSpecies=split /\_/, $tmp[2]; ## Need to modify if seperated by not underscore;
             	my $index=EBALib::CommonSubs::indexArray($tmp[2],@SpsArray);
	     	my $newIndex=$index+4; my $tarCoordinate="$tmp[5]--$tmp[6]";   ## Index is increase by 4 becuase fist three are species detail
			if (exists $doneHash{$newIndex}) { 
				my $oldTag=$doneHashPseudo{$newIndex};
				my @vals = values %doneHashPseudo;
				my @vals2=EBALib::CommonSubs::sortUniqueHash(@vals);
				## conditions : all Breaks and input file value is breaks or PsuedoBreaks and Break but infile values is Breaks
				if ((((scalar(@vals2) == 1) and ($vals2[0] eq "Break")) and ($tmp[8] eq "Break")) or (($tmp[8] eq "Break") and (scalar(@vals2) >=1))) {$multipleHit="YES";}
				my $oldCoordinate=$doneHash{$newIndex};
				my $newtarCoordinate = join (",", $oldCoordinate,$tarCoordinate);
				$doneHash{$newIndex}=$newtarCoordinate; 
				$doneHashPseudo{$newIndex}=$tmp[8]; 
				$tarCoordinate=$newtarCoordinate;
				}
			else 	{ $doneHash{$newIndex}=$tarCoordinate; $doneHashPseudo{$newIndex}=$tmp[8];}

		if ($tmp[8] eq "Break") { push (@allCoordinates, $tarCoordinate);}
		$matrix[$newIndex]=$tarCoordinate;
		$matrix[0]=$tmp[1];
		$matrix[1]=$tmp[0];
		$matrix[2]="$tmp[3]--$tmp[4]";
		}
my $decision=gapDecision(@allCoordinates);  ## check for the both gap conditions here !!!??????
if (($multipleHit eq "YES") or ($decision==0) or ($brkInfo eq "PseudoBreak")) { $matrix[3]="Gap"; } else { $matrix[3]="Breakpoints";}
undef @allCoordinates;   undef %doneHash;
return @matrix; 
}

##----------------------------------------------------------------------------------
        
sub gapDecision {  
my @allCoordinates=@_;  
my (@all_data, @val1,@val2); 
my $j1=0;my $j2=0;my $decision;

	foreach my $coordinate (@allCoordinates) {
		chomp ($coordinate);
		$coordinate=EBALib::CommonSubs::trim($coordinate);
		next if ($coordinate =~ m/^[0]$/);
		my @allCoordi=split /\,/, $coordinate;
		foreach my $cor (@allCoordi) {
			push (@all_data,$cor);
		}
	}

#	@all_data=clean_array(@all_data);

	if ($#all_data>1) {
	
		foreach my $f(@all_data) {
			next if (!$f);
			@val1=split /\-\-/, $f;

			foreach my $c(@all_data) {
				next if (!$c);
				@val2=split /\-\-/, $c;
				#print OUTFILE "$c\n";
				next if (($val1[0]==$val2[0]) && ($val1[1]==$val2[1]));
				my $OverRes = EBALib::CommonSubs::checkCorOverlaps ($val1[0],$val1[1],$val2[0],$val2[1]);
				if ($OverRes) { $j1++; } else 	{ $j2++; }
				}
		
			}
		}
	else	{ return $decision=1; }
		  	
if ($j2==0) {
return $decision=1;  #Breakpoints
}
else {
return $decision=0;  #Gaps
}
undef @all_data;

}



1;

