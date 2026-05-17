# $Id$ Find Reuse
# Perl module for EBA EBALib::FindReuse;
# Author: Jitendra Narayan <jnarayan81@gmail.com>
# Copyright (c) 2015 by Jitendra. All rights reserved.
# You may distribute this module under the same terms as Perl itself

##-------------------------------------------------------------------------##
## POD documentation - main docs before the code
##-------------------------------------------------------------------------##

=head1 NAME

EBALib::FindReuse  - DESCRIPTION of Object

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

package EBALib::FindReuse;
use strict;
#use warnings;
use File::Copy;
use List::Compare;
#use Term::ANSIColor;
use List::Util 'sum';
use List::Util qw(min max);

use Exporter;

our @EXPORT_OK = "generateReuse";

sub generateReuse {
my ($path, $spsNum, $threshold, $res, $engrave, $refName) = @_;
my $dir="$path/EBA_OutFiles";
EBALib::Messages::reuseMSG();

		 opendir(DIR, $dir) or die $!;
   		 while (my $file = readdir(DIR)) {
        		# We only want files
        		next unless (-f "$dir/$file");
        		# Use a regular expression to find files ending in .txt
        		next unless ($file =~ m/\.eba7$/); 
			my @new_file=split(/_/, $file);
				if ($new_file[0] ne "") { 
					my $InPutFile="$path/EBA_OutFiles/final_classify.eba7";
					my $OutPutFile= "$path/EBA_OutFiles/final_classify_reuse.eba8";
					createReuseFile($InPutFile, $OutPutFile, $spsNum, $threshold, $path, $res, $engrave, $refName);  ## We use $threshold=20 in bird dataset

				}				
		  }
closedir(DIR);
}

sub createReuseFile {

my ($InFile, $OutFile, $spsNum, $cutoff, $path, $res, $engrave, $refName)=@_;

open INFILE,  $InFile or die EBALib::Messages::failOp($InFile);
open OUTFILE, ">" , $OutFile or die EBALib::Messages::failOp($OutFile);

#print "$InFile\t$OutFile\t$spsNum\t$cutoff\n";

my %allClassificationHash; my @allSpeciesName;
open INCLASSFILE, 'classification.eba' or die EBALib::Messages::failOp("classification.eba");
while (<INCLASSFILE>) { chomp;  my $line= lc($_); $line=EBALib::CommonSubs::trim($line); #print "$line\n"; 
if (index($line,"#") == 0) { next; } # Lines starting with a hash mark are comments
if ($line =~ /^\s*$/) { next; }
my @tmp = split /\=/, lc($line); s{^\s+|\s+$}{}g foreach @tmp; push @allSpeciesName, $tmp[1]; $allClassificationHash{$tmp[0]} = $tmp[1]; } ## Store the classification file in Hash;

print OUTFILE "REName\tChrName\tClassName\tBrkCorStart\tBrkCorEnd\tDecision\tScores\tRatio\tSpeciesNumber\tGapNum\tBrkNum\tPercentageUsed\tRatioNew\n";  ## Now changed

$|++;
$/ = "\n";

while (<INFILE>) {
	my $line=lc($_); chomp $line; ## lowercase
	if ($line =~ /^\s*#/) { next; }
	next if $.== 1;   ### To remove the header ... if header is not there then we need to remove it
	my @tmp=split /\t/, $line; 
	my $end=scalar(@tmp);
	my $start=$spsNum+1;
	my @result = @tmp[$start..$end]; ## start from the chromosome column
	my @coordinates= split /\<\-\>/, $result[3];
	my @allClassGroup= split /\s+/, $result[5]; ## Need to modify if we add other seperators, it contains all the classification group names.
	my %allClassGroupHash;
	my $percentageUsed=(($result[10]-$result[11])/$result[10])*100;

	foreach (@allClassGroup) { my @GroupValues=split /\:/, $_;  $allClassGroupHash{$GroupValues[0]}=$GroupValues[1]; } ## store all the classification name and its scores.
	
		my @classNamePahala = split /\:/, $result[6]; ## first best classification
		my @classNameDusara = split /\:/, $result[7]; ## second best classification
		
		my @classValuesPahala; my @classValuesDusara;
		my $classValuePahala=$allClassificationHash{$classNamePahala[0]}; ## extracted first the class->species names.
		my $classValueDusara=$allClassificationHash{$classNameDusara[0]}; ## extracted the second class->species names.

		if (defined $classValuePahala and length $classValuePahala) { @classValuesPahala=split /\,/, $allClassificationHash{$classNamePahala[0]}; } else { push @classValuesPahala, $classNamePahala[0]; }
	
		if (defined $classValueDusara and length $classValueDusara) { @classValuesDusara=split /\,/, $allClassificationHash{$classNameDusara[0]}; } else { push @classValuesDusara, $classNameDusara[0]; }

		my $ListCompare = List::Compare->new(\@classValuesPahala, \@classValuesDusara);

    		my @intersection = $ListCompare->get_intersection; ## common between two.
    		my @union = $ListCompare->get_union;
		
		if($result[8] eq 'na' or $result[8] eq "") { $result[8]=0; EBALib::Messages::NA($.); } 

		if ((($result[8] <= $cutoff) and ($result[8] > 1)) and (scalar(@intersection) == 0)) { 

			checkReuse(\@coordinates, \@classNamePahala, \@classNameDusara, \@result, \@allClassGroup, \%allClassGroupHash, \%allClassificationHash, $percentageUsed, $cutoff, 1, $refName); ## flag 1 if non-overlapping means REUSE
			
			}
		elsif ($result[8] == 1) {
			print OUTFILE "$refName\t$result[0]\t$classNamePahala[0]:$classNameDusara[0]\t$coordinates[0]\t$coordinates[1]\tUncertain\t$classNamePahala[1]:$classNameDusara[1]\t$result[8]\t$result[10]\t$result[11]\t$result[12]\t$percentageUsed\t$result[8]\n";
			}
		else {
			my $signal=checkReuse(\@coordinates, \@classNamePahala, \@classNameDusara, \@result, \@allClassGroup, \%allClassGroupHash, \%allClassificationHash, $percentageUsed, $cutoff, 0, $refName); ## Overlapping
			if ($signal == 0) {			
				print OUTFILE "$refName\t$result[0]\t$classNamePahala[0]\t$coordinates[0]\t$coordinates[1]\tUnique\t$classNamePahala[1]\t$result[8]\t$result[10]\t$result[11]\t$result[12]\t$percentageUsed\t$result[8]\n";
				}
 
			}
 
	undef @classValuesDusara; 
	#}
#undef @allClassGroup;
}
close INFILE or die EBALib::Messages::failCl("$InFile");
close OUTFILE or die EBALib::Messages::failCl("$OutFile");

## Generate a result files in ResultFiles folder
generateFinal("$path/EBA_OutFiles/final_classify.eba7", "$path/ResultFiles/Result_$res.final", 1, $spsNum, $engrave); 
generateFinal("$path/EBA_OutFiles/final_classify_reuse.eba8", "$path/ResultFiles/ResultReuse_$res.final", 2, $spsNum, $engrave);

#Create a final STATS
FinalStater("$path/EBA_OutFiles/final_classify.eba7", "$path/ResultFiles/Result_$res.stats", $spsNum); 

# Move the files
#copy ("$path/EBA_OutFiles/final_classify.eba7", "$path/ResultFiles/") or die "Failed to copy $_: $!\n";
#copy ("$path/EBA_OutFiles/ResultReuse.final", "$path/ResultFiles/") or die "Failed to copy $_: $!\n";
}


#Create a final Stats of EBRs
sub FinalStater{
my ($filename,$outfile,$spsNum)=@_;
open(FFILE,"<$filename") or die "Can't open configuration file $filename.";
my %countSTAT;
my %countSTAT1;
my @ff=<FFILE>;
shift @ff; # Remove the first element HEADER and throw it away
foreach my $line (@ff) {
	my @ln=split("\t",$line);
	my $finalCor=$spsNum+7; # 7 the column to extends!!
	my @cla =split('\:', $ln[$finalCor]);
	my $fiCor=$spsNum+9;
	if ($ln[$fiCor] == 1) { $countSTAT1{$cla[0]}++; next;}
	$countSTAT{$cla[0]}++;
}
close CFGFILE;
open(OUTFILEWA, ">", $outfile) or die EBALib::Messages::failOp($outfile); 
if (-z $outfile) {print OUTFILEWA "ClassificationGroup\tEBRs\n";}
foreach my $str (sort keys %countSTAT) { print OUTFILEWA "$str\t$countSTAT{$str}\n"; }
#Print the EBRs with 1s
#foreach my $str1 (sort keys %countSTAT1) { print OUTFILEWA "$str1\t$countSTAT1{$str1}\n"; }
close OUTFILEWA;
}

#Check the re-use cases
sub checkReuse {

my ($coordinates_ref,$classNamePahala_ref, $classNameDusara_ref, $result_ref, $allClassGroup_ref, $allClassGroupHash_ref, $allClassificationHash_ref, $percentageUsed, $cutoff, $flag, $refName)=@_;
my @coordinates=@$coordinates_ref; 
my @classNamePahala=@$classNamePahala_ref; 
my @result=@$result_ref; 
my @allClassGroup=@$allClassGroup_ref; 
my %allClassificationHash=%$allClassificationHash_ref; 
my @classNameDusara=@$classNameDusara_ref;
my %allClassGroupHash=%$allClassGroupHash_ref;

my @storeName; 
my @storeScore; 
my @storeRatio; 
my @nextScore; 
my $signal=0;
#my $nextMaximum=0;

	if($flag == 1) { ## if it reuse		
		push @storeName, $classNamePahala[0]; push @storeScore, $classNamePahala[1]; push @storeRatio, $result[8]; ## Store only once becuase we need to print once
		push @storeName, $classNameDusara[0]; push @storeScore, $classNameDusara[1]; 
		my @allGroupNameArray; foreach my $keys (keys %allClassificationHash) { push(@allGroupNameArray,$keys);} ## stored classification.eba keys for future use

		foreach my $key (sort { $allClassGroupHash{$a} <=> $allClassGroupHash{$b} || $a cmp $b } (keys %allClassGroupHash) ) {
				
			next if (($classNamePahala[1] == $allClassGroupHash{$key}) or ($classNameDusara[1] == $allClassGroupHash{$key}));  ## Becuase we have allready store them if non-overlapping[REUSE]... 
					
			my @cNP = split /\:/, $result[6]; ## this is the .eba7 breaks best score column
		
			my @cVP; my @cVD; 
			my $cVP=$allClassificationHash{$cNP[0]}; 
			my $cVD=$allClassificationHash{$key}; ## The second 

			if (defined $cVP and length $cVP) { @cVP=split /\,/, $allClassificationHash{$cNP[0]}; } else { push @cVP, $cNP[0]; }
	
			if (defined $cVD and length $cVD) { @cVD=split /\,/, $allClassificationHash{$key}; } else { push @cVD, $key; }

			my $LC = List::Compare->new(\@cVP, \@cVD);

    			my @i = $LC->get_intersection;
    			my @u = $LC->get_union;

			next if $allClassGroupHash{$key} == 0;
			my $r=$cNP[1]/$allClassGroupHash{$key};
			my $presentIn=EBALib::CommonSubs::isInList($key, @allGroupNameArray);
			if ((($r <= $cutoff) and ($r > 1)) and (scalar(@i) == 0) and ($presentIn == 0) ) { 

				push @storeName, $key; push @storeScore, $allClassGroupHash{$key}; push @storeRatio, $r;
				$signal=1;
				}
			else { push @nextScore, EBALib::CommonSubs::expand($allClassGroupHash{$key});  }


		undef @cVD; undef @cVP; 
		}

		my $shiftedName=""; 
		my $shiftedScore=""; 
		my $shiftedRatio=""; 		
		
		for my $num (0 .. $#storeName)  { 
			
			if ($num > 0) { push @storeName, $shiftedName; push @storeScore, $shiftedScore; push @storeRatio, $shiftedRatio; } ## push the first name at the last of array ..
			my $secondIsPresent=EBALib::CommonSubs::isInList($classNameDusara[0], @allGroupNameArray);

			if (($num == 1) and ($secondIsPresent == 1) and (scalar(@storeName) > 2)) { ## It delete the second best score in second round of printing [ if reuse is more than 2] in outfile.. 
				delete $storeName[0]; delete $storeScore[0]; delete $storeRatio[0]; @storeName = grep{$_} @storeName;  @storeScore = grep{$_} @storeScore;  @storeRatio = grep{$_} @storeRatio;
				}
				
			my $myClassNameMain= join(":",@storeName);
			my $myClassScoreMain= join(":",@storeScore);
			my $myClassRatioMain= join(":",@storeRatio);
			my $sum=sum(@storeScore);
			my $sumRatio=sum(@storeRatio);
			my $averageScore=$sum/scalar(@storeScore);
			my $nextMaximum=max(@nextScore);
			my $minScore=min(@storeScore);
			my $diffValue=$minScore-$nextMaximum;
			my $newScoreReuse;
			if ($diffValue < 0) { $newScoreReuse=$sumRatio/scalar(@storeRatio);} 
			else { if (!$nextMaximum) {$newScoreReuse = "Undef";} else {$newScoreReuse=$averageScore/$nextMaximum; } } # Sometime the nextMaximum==0 !!

			print OUTFILE "$refName\t$result[0]\t$myClassNameMain\t$coordinates[0]\t$coordinates[1]\tReuse\t$myClassScoreMain\t$myClassRatioMain\t$result[10]\t$result[11]\t$result[12]\t$percentageUsed\t$newScoreReuse\n"; 

			$shiftedName=shift @storeName;
			$shiftedScore=shift @storeScore;
			$shiftedRatio=shift @storeRatio;
			
		}
	undef @storeName; undef @storeScore; undef @nextScore; undef $shiftedName; undef $shiftedScore; undef $shiftedRatio; 
	} ## if ends here

##-----------------------------------

	elsif($flag == 0) {

		my @allGroupNameArray; foreach my $keys (keys %allClassificationHash) { push(@allGroupNameArray,$keys);} ## stored classification.eba keys for future use

		foreach my $key (sort { $allClassGroupHash{$a} <=> $allClassGroupHash{$b} || $a cmp $b } (keys %allClassGroupHash) ) {
				
			next if (($classNamePahala[1] == $allClassGroupHash{$key}) or ($classNameDusara[1] == $allClassGroupHash{$key}));  ## Becuase we have allready store them if non-overlapping[REUSE]... 
					
			my @cNP = split /\:/, $result[6]; ## this is the .eba7 breaks best score column
		
			my @cVP; my @cVD; 
			my $cVP=$allClassificationHash{$cNP[0]}; 
			my $cVD=$allClassificationHash{$key}; ## The second 

			if (defined $cVP and length $cVP) { @cVP=split /\,/, $allClassificationHash{$cNP[0]}; } else { push @cVP, $cNP[0]; }
	
			if (defined $cVD and length $cVD) { @cVD=split /\,/, $allClassificationHash{$key}; } else { push @cVD, $key; }

			my $LC = List::Compare->new(\@cVP, \@cVD);

    			my @i = $LC->get_intersection;
    			my @u = $LC->get_union;

			next if $allClassGroupHash{$key} == 0;
			my $r=$cNP[1]/$allClassGroupHash{$key};
			my $presentIn=EBALib::CommonSubs::isInList($key, @allGroupNameArray);
			if ((($r <= $cutoff) and ($r > 1)) and (scalar(@i) == 0) and ($presentIn == 0) ) { 

				push @storeName, $key; push @storeScore, $allClassGroupHash{$key}; push @storeRatio, $r;
				$signal=1;
				}
			else { push @nextScore, EBALib::CommonSubs::expand($allClassGroupHash{$key});  }


		undef @cVD; undef @cVP; 
		if($classNameDusara[1] == $allClassGroupHash{$key}) { unshift @nextScore, $allClassGroupHash{$key};}  ## Becuase first and second is not reuse in this case .. we need to store second score
		}
		
		if($signal == 1) { unshift @storeName, $classNamePahala[0]; unshift @storeScore, $classNamePahala[1]; unshift @storeRatio, $result[8]; } ## This work only when reuse in first unique condition are checked .. !!! I store the pahala class info only when there is some reuse ... by looking at the signal
		my $shiftedName=""; 
		my $shiftedScore=""; 
		my $shiftedRatio=""; 		

		for my $num (0 .. $#storeName)  { 

			if ($num > 0) { push @storeName, $shiftedName; push @storeScore, $shiftedScore; push @storeRatio, $shiftedRatio; } ## push the first name at the last of array ..
				
			my $myClassNameMain= join(":",@storeName);
			my $myClassScoreMain= join(":",@storeScore);
			my $myClassRatioMain= join(":",@storeRatio);
			my $sum=sum(@storeScore);
			my $sumRatio=sum(@storeRatio);
			my $averageScore=$sum/scalar(@storeScore);
			my $nextMaximum=max(@nextScore);
			my $minScore=min(@storeScore);
			my $diffValue=$minScore-$nextMaximum;
			my $newScoreReuse;
			if ($diffValue < 0) { $newScoreReuse=$sumRatio/scalar(@storeRatio);}
			else { if (!$nextMaximum) {$newScoreReuse = "Undef";} else {$newScoreReuse=$averageScore/$nextMaximum; } } # Sometime the nextMaximum==0 !!

			print OUTFILE "$refName\t$result[0]\t$myClassNameMain\t$coordinates[0]\t$coordinates[1]\tReuse\t$myClassScoreMain\t$myClassRatioMain\t$result[10]\t$result[11]\t$result[12]\t$percentageUsed\t$newScoreReuse\n"; 

			$shiftedName=shift @storeName;
			$shiftedScore=shift @storeScore;
			$shiftedRatio=shift @storeRatio;
			
		}
	undef @storeName; undef @storeScore; undef @nextScore; undef $shiftedName; undef $shiftedScore; undef $shiftedRatio; #undef @hogaya;
	return $signal;	
		
	}
	
} ## Reuse function ends here
	
sub generateFinal {
use Math::Round;
    my ($file1, $file2, $fnum, $spsNum, $engrave) = @_;
    my @lines; my @tmpLine;
    open(FILE1, $file1) or die EBALib::Messages::failOp($file1);
    open(FILE2, ">", $file2) or die EBALib::Messages::failOp($file2);
    while(<FILE1>){
        chomp; #remove newline
        s/(^\s+|\s+$)//g; # strip lead/trail whitespace
        next if /^$/;  # skip blanks
	@tmpLine = split(/\t/, $_);
        #push @lines, $_; 
	if (($fnum==1) and ($engrave == 0)){
		for (my $aa=0; $aa<=($spsNum-1); $aa++) { 
			#print FILE2 "$tmpLine[$aa]\t"; 
			}
		if ($. == 1) { print FILE2 "Chromosome\tWidest EBR interval start (bp)\tWidest EBR interval end (bp)\tNarrowest  EBR interval start (bp)\tNarrowest  EBR interval start (bp)\tSpecies containing the EBR\tClassification group:assignment score\tHighest probability classification\tSecond highest probability classification\tRatio between the first and second probabilities\tNo. species in the dataset\tNo. species with gap\tNo. species with the EBR\tPercentage informative species\tError probability\n";  next;}
		my @narrowCoordi=split(/\<\-\>/, $tmpLine[$spsNum+4]);
		my $perUsed=(($tmpLine[$spsNum+11]-$tmpLine[$spsNum+12])/$tmpLine[$spsNum+11])*100;
		my $roundRatio;
		if(($tmpLine[$spsNum+9]) ne "NA") { $roundRatio=sprintf("%.2f", $tmpLine[$spsNum+9]);} else { $roundRatio="NA";}
		my $newperUsed=round($perUsed);
		my $pval;
		if ($roundRatio eq "NA" ) { $roundRatio="Undef";} else { $pval=1/$roundRatio;}
		my @newWide=split (/\-\-/, $tmpLine[$spsNum+2]);
		print FILE2 "$tmpLine[$spsNum+1]\t$newWide[0]\t$newWide[1]\t$narrowCoordi[0]\t$narrowCoordi[1]\t$tmpLine[$spsNum+5]\t$tmpLine[$spsNum+6]\t$tmpLine[$spsNum+7]\t$tmpLine[$spsNum+8]\t$roundRatio\t$tmpLine[$spsNum+11]\t$tmpLine[$spsNum+12]\t$tmpLine[$spsNum+13]\t$newperUsed\t$pval\n";
		
	}
	elsif (($fnum==1) and ($engrave == 1)){
		for (my $aa=0; $aa<=($spsNum-1); $aa++) { print FILE2 "$tmpLine[$aa]\t"; }
		if ($. == 1) { print FILE2 "Chromosome\tWidest EBR interval start (bp)\tWidest EBR interval end (bp)\tNarrowest  EBR interval start (bp)\tNarrowest  EBR interval start (bp)\tSpecies containing the EBR\tClassification group:assignment score\tHighest probability classification\tSecond highest probability classification\tRatio between the first and second probabilities\tNo. species in the dataset\tNo. species with gap\tNo. species with the EBR\tPercentage informative species\tError probability\n";  next;}
		my @narrowCoordi=split(/\<\-\>/, $tmpLine[$spsNum+4]);
		my $perUsed=(($tmpLine[$spsNum+11]-$tmpLine[$spsNum+12])/$tmpLine[$spsNum+11])*100;
		my $roundRatio;
		if(($tmpLine[$spsNum+9]) ne "NA") { $roundRatio=sprintf("%.2f", $tmpLine[$spsNum+9]);} else { $roundRatio="NA";}
		my $newperUsed=round($perUsed);
		my $pval;
		if ($roundRatio eq "NA" ) { $roundRatio="Undef";} else { $pval=1/$roundRatio;}
		my @newWide=split (/\-\-/, $tmpLine[$spsNum+2]);
		print FILE2 "$tmpLine[$spsNum+1]\t$newWide[0]\t$newWide[1]\t$narrowCoordi[0]\t$narrowCoordi[1]\t$tmpLine[$spsNum+5]\t$tmpLine[$spsNum+6]\t$tmpLine[$spsNum+7]\t$tmpLine[$spsNum+8]\t$roundRatio\t$tmpLine[$spsNum+11]\t$tmpLine[$spsNum+12]\t$tmpLine[$spsNum+13]\t$newperUsed\t$pval\n";
		
	}
	elsif($fnum==2){
		if ($. == 1) { print FILE2 "Reference genomme name\tReference Chr.\tFinal classification\tNarrowest  EBR interval start (bp)\tNarrowest  EBR interval end (bp)\tType\tProbabilities\tRatio between the first and second classification probabilities\tNo. species in the dataset\tNo. species with gap\tNo. species with the EBR\tPercentage informative species\tRatio between the first and second classification probabilities adjusted for reuse EBRs\n"; next;}
		my $newroundRatio;
		if(($tmpLine[12]) ne "NA") { $newroundRatio=sprintf("%.2f", $tmpLine[12]);} else { $newroundRatio="NA";}
		my $newperUsed2=round($tmpLine[11]);
		print FILE2 "$tmpLine[0]\t$tmpLine[1]\t$tmpLine[2]\t$tmpLine[3]\t$tmpLine[4]\t$tmpLine[5]\t$tmpLine[6]\t$tmpLine[7]\t$tmpLine[8]\t$tmpLine[9]\t$tmpLine[10]\t$newperUsed2\t$newroundRatio\n";
	}
    }
    close FILE1;
    close FILE2;
}

1;
