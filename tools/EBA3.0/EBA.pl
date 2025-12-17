#!/usr/bin/env perl
#BEGIN { $^W = 1 }

# $Id: EBA.pl 01 2015-06-15 23:25:00Z Jitendra $

##---------------------------------------------------------------------------##
##  File: EBA.pl
##
##  Author:
##        Jitendra <jnarayan81@gmail.com>, Denis <dmlarkin@gmail.com>
##
##  Description:
##
#******************************************************************************
#* Copyright (C) 2015 Jitendra Narayan <jnarayan81@gmail.com>
#* This work is distributed under the license of Perl itself.
###############################################################################

##---------------------------------------------------------------------------##
## Module dependencies
##---------------------------------------------------------------------------##

#!/usr/bin/perl

use strict;
use warnings;

use FindBin;
use lib "$FindBin::Bin";
use English;
use FileHandle;
use Getopt::Long;
use Pod::Usage;
use File::Path;
use Cwd;
use File::Find;
use File::Copy;

use EBALib::StoreSpecies;
use EBALib::CheckData;
use EBALib::CalculateBeta;
use EBALib::ModifyBeta;
use EBALib::BreaksFinder;
use EBALib::BreaksAmongstSpecies;
use EBALib::ConCatFile;
use EBALib::BreaksMatrix;
use EBALib::BreaksScoring;
use EBALib::BreaksScoring2;
use EBALib::EnterScore;
use EBALib::ConCatAll;
use EBALib::CreateFinal;
use EBALib::PoissonMethod;
use EBALib::ClassifyBreakpoints;
use EBALib::CheckClassification;
use EBALib::ModifyClassification;
use EBALib::MergeResolution;
use EBALib::FindReuse;
use EBALib::Messages;
use EBALib::CommonSubs;
use EBALib::Draw::drawPieChart;
use EBALib::Draw::drawPieChartFinal;
use EBALib::Draw::drawPieChartBreaksFinal;
use EBALib::Draw::drawCumulatedBar;
use EBALib::Draw::drawCumulatedStackedBar;
use EBALib::Draw::drawCumulatedStackedBarMerged;
use EBALib::Draw::drawPieChartBreaksFinalMerged;
use EBALib::Draw::drawBrkChrLineGraphFinal;
use EBALib::FindReuseMerge;
use EBALib::Draw::drawPieChartFinalMerged;
use EBALib::Draw::drawBrkChrLineGraphFinal2;
use EBALib::Draw::drawBreakpointChrGraphFinal;
use EBALib::Draw::drawBreakpointChrGraphFinal2;
use EBALib::Draw::drawCumulatedStackedBarAll;
use EBALib::Safai;

use File::Basename qw(dirname);
use Cwd qw(abs_path);
use FindBin;
use lib "$FindBin::RealBin/..";
use lib qw(/home/bajay/Ajay_Bhatia/home_scripts/EBA_folder/eba/EBA3.0/bin);
use lib '/home/bajay/Ajay_Bhatia/lab1/Rotifer_study2/tools/EBA3.0/bin/EBALib';


my $myCurDir = getcwd;
my $VERSION = 'v3.O';
print <<'WELCOME'; 
	  ______ _____     ___   
	 | ____| |  _ \	   / \  
	 |  _|   | |_) |  / _ \   	
	 | |___  | |_) / / ___ \  
	 |_____| |____/ /_/   \ \  
	 EBA v3.0 [Sept 13, 2016]  - [+]

Citation - automated definition and classification of Evolutionary Breakpoints Regions
Paper URL

License: Creative Commons Licence
Bug-reports and requests to: jnarayanATgmail.com and dmlarkinATgmail.com
WELCOME

EBALib::CommonSubs::returnOS($VERSION);
$|++;

# Quick help
unless (@ARGV) { 
	# when no command line options are present
	# print SYNOPSIS
	pod2usage( {
		'-verbose' => 0, 
		'-exitval' => 1,
	} );
}

my (
	$number,
	$classify, 
	$reference, 
	$dir,
	$outdir, 
	$force, 
	$beta, 
	$validate, 
	$lineage, 
	$prime,
	$merge, 
	$exclude,
	$scrutiny,
	$increase,
	$threshold,
	$verbose,
	$logfile,
	$keep,
	$engrave,
);

# Default option setting for EBA tool
$verbose=0; 		# Verbose set to 0;
$increase=0; 		# The default breakpoint size increase is set to zero;
my %options = ();

GetOptions(
	\%options,
    	'number|n=i'    => \$number,        	## Number of species on are working with
    	'dir|d=s' 	=> \$dir,           	## In directory name
	'outdir|o=s' 	=> \$outdir,           	## Out Directory name
    	'beta|b=s' 	=> \$beta,          	## Beta Calculation "yes" or "no"
    	'force|f' 	=> \$force,         	## if "yes" work on unlimited species and resolutions
    	'classify|c=s' 	=> \$classify,  	## "yes" or "no" ## classify the species or not
    	'validate|v' 	=> \$validate,  	## "yes" or "no" ## validate the species or not
	'lineage|l' 	=> \$lineage,  		## "yes" or "no" ## print lineage or not
	'prime|p=n'	=> \$prime,  		## primary resolution name from you available list ## e.g 100 
	'merge|m'	=> \$merge,  		## Merging all resolutions using PRIME resolution as a base
	'increase|i=n'	=> \$increase,  	## the increment size provided in base pair ## e.g 5000
	'threshold|t=n'	=> \$threshold,  	## the threshold to filter the reuse breakpoints.
	'exclude|x'	=> \$exclude,  		## Exclude the classification group if only one species is there to decide the order/group.
	'keep|k' 	=> \$keep,		## If use want to keep the intermediate file provide "yes" ... default will delete all the intermediate files.
	'engrave|e' 	=> \$engrave,		## print the breaks detail of all studied species.
	'scrutiny|s' 	=> \$scrutiny,		## scrutiny the classification file
    	'reference|r=s' => \$reference, 	## Name of your reference species 
    	'help|?|h!'     => sub { EBALib::Messages::EBAWelcome($VERSION) },
   	'who|w!'     	=> sub { EBALib::Messages::EBAWho($VERSION) },
	'verbose' 	=> \$verbose,
    	'logfile=s' 	=> \$logfile,		## The name of the log file .. We will add .log extension in it.
	
) or die EBALib::Messages::ManualHelp();

## make sure everything passed was peachy
#EBALib::Messages::checkParameters(\%options);

#Store in current dir
if (!$outdir) { $outdir = "EBA_OUT";}
EBALib::CommonSubs::dircopy($dir,$outdir);

#Check the taxdump folder, and if not there download it
#use LWP::Simple;
#my $url = 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz';
#my $taxdump = 'taxdump.tar.gz';
#if (!$taxdump){ getstore($url, $taxdump);}
#use Archive::Extract;
### build an Archive::Extract object ###
#my $ae = Archive::Extract->new( archive => 'taxdump.tar.gz' );
### extract to cwd() ###
#my $ok = $ae->extract or die $ae->error;

if ((!$number) or (!$dir) or (!$reference) or (!$prime) or (!$threshold)) { EBALib::Messages::printUsage(); }

# If the beta score is not calculated then user need to provide their own beta score in an acceptable format.
if ($beta) { 
	EBALib::Messages::BetaUpdate($myCurDir); 
	if (-e $beta) { 
		EBALib::Messages::BetaExists() } 
	else { 
		EBALib::Messages::BetaMissed($myCurDir)
	}
}

# If the classification is not calculated then user need to provide their own classification file in acceptable format.
if ($classify) { 
	EBALib::Messages::ClassifyNoNCBI($myCurDir);
	if (-e $classify) { 
		EBALib::Messages::ClassifyCustom() } 
	else { 
		EBALib::Messages::ClassifyProvide($myCurDir);
	}
}

# Currently EBA restricted to only TEN species. However, for internal user they can use -f or -force=yes option to forcefully run it.
if (($number > 10) and (!$force)) { EBALib::Messages::EBAMax ($VERSION);}
else { 
	my $startTime=time; 	# Set the starting time, use for calculating the time taken in analysis.
	if (defined $logfile) {
  		$verbose = 1;
		my $newlogfile="$logfile.log"; 	# The default logfile name if call at command prompt.
  		open my $log_handle, '>', $newlogfile or die "Could not open $newlogfile";
  		# select makes print point to LOGFILE 
  		select($log_handle);
	}	
	
	if(($number and $dir and $reference and $prime and $threshold) eq "") { print "ERROR : Wrong Usage "; EBALib::Messages::printUsage();}
	## Delete merge folder if exist
	rmtree(["$outdir/Merge"]); 
	
	## Read the resolutions folders, files and validate it for the EBA.
	if (!$validate) { EBALib::CheckData::verifyData($outdir, $number); } else { $validate=1;}

	my $parent = "./$outdir";  #The main folder which contains all the resolution folders

	my ($par_dir, $sub_dir);

	opendir ($par_dir, $parent) or die "cannot open $parent directory";
	while (my $sub_folders = readdir($par_dir)) {
    	next if ($sub_folders =~ /^\.\.?$/);  # skip . and ..
    	my $path = $parent . '/' . $sub_folders;
    	next unless (-d $path);   # skip anything that isn't a directory
		# my @name = split(/\//, $path);
		#my $currentDir = `pwd`; for linux ....
		my $currentDir = getcwd;
		print "\nWorking in $currentDir$path\n";
		my @res = split(/\//, $path);
		my $resName=$res[-1]; ## Resolution name
		my ($refName, $fileNames_ref, $allRes_ref)=calculate($path, $increase);
		my @fileNames=@$fileNames_ref;
		my @allRes=@$allRes_ref;

		if (!$classify) {
			EBALib::Messages::WorkOnClass();
			EBALib::ClassifyBreakpoints::classifyEBA($reference); 
			$classify=1;
		}
		
		if (!$scrutiny) { 
			EBALib::Messages::ParseClass();
			EBALib::CheckClassification::verifyClassification('classification.eba', $reference); $scrutiny=1;}
		
		if ($exclude) { 
			EBALib::Messages::ExludeClass();
			EBALib::ModifyClassification::alterClassification('classification.eba'); $exclude=0;}

		EBALib::BreaksAmongstSpecies::breakpointsAmongstSpecies($path);
		EBALib::BreaksMatrix::breakpointsTable($path);
		EBALib::Draw::drawPieChart::drawPie($path, $resName, $number);
		EBALib::Draw::drawCumulatedBar::drawCumBar($path, $resName, $number);
		EBALib::BreaksScoring::breakpointScoring1($path,@fileNames);
		EBALib::BreaksScoring2::breakpointScoring2($path);
		EBALib::EnterScore::enterScore($path,@fileNames);
		EBALib::ConCatAll::concatAll($path);
		EBALib::CreateFinal::generateFinal($path); 
		
		if (!$beta) { 
			EBALib::Messages::CalBS();
			EBALib::CalculateBeta::betaCal($outdir, $increase, $number, $resName);

			use strict;
			use File::Find;
			find(\&dir_names, "$outdir");
			my (@allDIR, @finalDIR);
			sub dir_names {
 				@allDIR="$File::Find::dir" if(-f $File::Find::dir,'/');
				foreach(@allDIR) { 
					my @resDIR = split(/\//, $_);
					if ($resDIR[-1] =~ /^[+-]?\d+$/ ) { ## If number
					#print $resDIR[-1];
					push @finalDIR, $resDIR[-1];
					}
				}
			}

			@finalDIR=EBALib::CommonSubs::uniq(@finalDIR);
			my @finalDIR_sorted = sort { $a <=> $b } @finalDIR;
			#foreach (@finalDIR_sorted) { print "$_\n";}
			my $resName=EBALib::CommonSubs::isInList($prime, @finalDIR_sorted);
			if ($resName == 0) { EBALib::Messages::ResRange($prime)}
			EBALib::ModifyBeta::alterBeta("betaScore", $number, $finalDIR_sorted[-1]);
			$beta=1; 
		}
		EBALib::PoissonMethod::generateFinalPoisson($path, $number, $lineage); ## Read number of species from command line 
		EBALib::Draw::drawPieChartFinal::drawPie($path, $number, $resName);
		EBALib::FindReuse::generateReuse($path, $number, $threshold, $resName, $engrave, $refName);
		EBALib::Draw::drawCumulatedStackedBar::drawStackBar($path,$resName, $number);
		EBALib::Draw::drawPieChartBreaksFinal::drawBreaksPie($path, $number,$resName);
		EBALib::Draw::drawBrkChrLineGraphFinal::breakpointGraphFinal($path,$resName, $number); ## To generate a line graph for merged final file
		EBALib::Draw::drawBrkChrLineGraphFinal2::breakpointGraphFinal($path,$resName, $number); ## To generate a line graph for merged final file
			
	# To print the time taken to accomplish the task
	my $tim=(time - $startTime) / 60;
	#printf "\nTotal time taken is $tim minutes\n";
	EBALib::Messages::TimeTaken($startTime);
	$startTime=time; ## Reset the time
	}
	closedir($par_dir);

# Merge stuff here
if (!$merge) {
   my $curDir = getcwd;
   EBALib::Messages::MergeMSG($curDir, $outdir);
   EBALib::MergeResolution::mergeAll($outdir, $number, $prime, $lineage); ## Merging all the resolutions
   EBALib::FindReuseMerge::generateReuseMerge($number, $threshold, $engrave); ## Calculate the merged resolutions
   EBALib::Draw::drawCumulatedStackedBarMerged::drawStackBarMerged($number); # Merged
   EBALib::Draw::drawPieChartFinalMerged::drawPieChartMerge($number);
   EBALib::Draw::drawPieChartBreaksFinalMerged::drawBreaksPieMerged($number);
   
   EBALib::Draw::drawBreakpointChrGraphFinal::breakpointGraphFinal("Result_Reuse_Merge.final", "Merge", $number); ## To generate a line graph for merged final file
   EBALib::Draw::drawBreakpointChrGraphFinal2::breakpointGraphFinal("Result_Reuse_Merge.final", "Merge", $number); ## To generate a line graph for merged final file
	unlink('Result_Merge2.final');
   $merge=1;
   }

# Make Merge and move all data.
mkdir ("$outdir/Merge", 0777) or print "$!\n";
moveAll("$outdir/Merge/", 'data'); 	##To move all *.data files
EBALib::Draw::drawCumulatedStackedBarAll::drawStackBarMergedAll($outdir, $number); ## It does not require to the same address.
my $curDir = getcwd;
deleteColumn ("$curDir/Result_Merge3.final", "$curDir/final_classify.final", $number);
moveAll("$outdir/Merge/", 'gif'); 		##To move all *.gif files
moveAll("$outdir/Merge/", 'final'); 	##To move all *.gif files
#print "The merged results and graph datasets are moved to $curDir/$outdir/Merge folder\n";
unlink glob "*.data"; unlink glob "*.tmp"; unlink glob "*.eba0";
## To clean up at the end of the program.
if (!$keep) { EBALib::Messages::cleanIn(); EBALib::Safai::cleanUp($outdir);} 

}

# Subroutines starts here -------------

sub calculate {
my ($path, $increase)=@_;
my @speciesNames; my $sub_dir; my $refName; my @allRes;
if(!$increase) { $increase=0;} ## the default increament for the breakpoint size is 0;

	opendir($sub_dir, $path); 
  	
	## Remove folder if already exist. 
	if (rmtree(["$path/EBA_OutFiles", "$path/EBA_ImageFiles", "$path/ResultFiles"],)) { 
		EBALib::Messages::deleteIn();
	} 
		#Create files in current folder		
		EBALib::Messages::createMSG();
		mkdir ("$path/EBA_OutFiles", 0777) or print "$!\n";
		mkdir ("$path/EBA_ImageFiles", 0777) or print "$!\n";
        	mkdir ("$path/ResultFiles", 0777) or print "$!\n";

    	while (my $file = readdir($sub_dir)) {
        	next unless $file =~ /\.txt?$/i;
        	my $full_path = $path . '/' . $file;
		my @first_name = split(/_/, $file); ## the input file should be separated with underscore
		#if ($first_name[0] eq /^\s*$/) { sensorium(); }    # Warn and make sense of errors.
		push (@speciesNames, $first_name[0]);
	
		$refName=storeRef($path, $file);
		my @res = split(/\//, $full_path);
		push (@allRes, $res[2]);
		EBALib::BreaksFinder::defineBreakpoints($path,$file,$increase);	 # The function is in Perl package BreaksFinder;
		EBALib::StoreSpecies::storeSpecies($path, $file);
	}
	EBALib::ConCatFile::concatNsps($path);
	closedir($sub_dir);
my @allRes_sorted = sort { $a <=> $b } @allRes;
return $refName, \@speciesNames, \@allRes_sorted;
}

sub storeRef {
my ($path, $file)=@_;
my $InFile="$path/$file";
my @refName;
open INFILE,  $InFile or die EBALib::Messages::failOp($InFile);
	while (<INFILE>) { chomp;  my $line= lc($_); $line=EBALib::CommonSubs::trim($line);
	my @tmp = split /\t/, lc($line); s{^\s+|\s+$}{}g foreach @tmp; push @refName, $tmp[0];} 
close INFILE or die EBALib::Messages::noClose();
return $refName[0];
}

#Move all files and delete
sub moveAll {
my ($copytopath, $ext)=@_;
$ENV{'Merge'}="$copytopath";
my @files = glob("*.$ext"); ## We can add path if required glob("$PATH1/*.data");
	for my $file (@files) {
    	copy("$file", $ENV{'Merge'}) or die EBALib::Messages::fail(); ## We can add the path here as well
	unlink ($file);
	}
}

sub deleteColumn {
use autodie;
my ($infile, $outfile, $spsNum) = @_;

open (my $info, '<', $infile) or die EBALib::Messages::failOp($infile);
open (my $out, '>', $outfile) or die EBALib::Messages::failCl($outfile);

my @remove;

while (<$info>) {
  chomp;
  my @data = split /\t/, $_;
  unless (@remove) {
    @remove = ($spsNum+1, $spsNum+2, $spsNum+3, $spsNum+6); ## Provide the indexnumber to delete
  }
  splice @data, $_, 1 for reverse @remove;
  print $out join("\t", @data), "\n";
}
close $info;
close $out;
unlink($infile);
}

### POD Documentation

__END__

=head1 NAME

EBA.pl	- Script to automated definition and classification of evolutionary breakpoints regions from a large number of genomes based on the species phylogenetic relationship.

=head1 SYNOPSIS

perl EBA.pl -n <number> -d <dir> -r <refname> -p <prime> -t <number> 

Mandatory parameters:

  -n <number> 	Provide the number of species that you are going to work with.

  -d <dir> 	The name of directory which contains the list of HSBs. 

  -r <refname> 	Name of your reference species [ Must be scirentific name ].

  -p <prime> 	Provide the primary resolution name [ the resolution name should be numeric ].

  -t <number> 	The threshold value for reuse breakpoint filtration. 

Try -h for more detail.

=head1 OPTIONS
