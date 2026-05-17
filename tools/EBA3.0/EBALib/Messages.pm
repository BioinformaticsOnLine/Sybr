# $Id$ All messages
# Perl module for EBA EBALib::Messages;
# Author: Jitendra Narayan <jnarayan81@gmail.com>
# Copyright (c) 2015 by Jitendra. All rights reserved.
# You may distribute this module under the same terms as Perl itself

##-------------------------------------------------------------------------##
## POD documentation - main docs before the code
##-------------------------------------------------------------------------##

=head1 NAME

EBALib::Messages  - DESCRIPTION of Object

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

package EBALib::Messages;
use warnings::register;

sub ManualHelp {
print("Unrecognized option(s)!! Please check manual OR Try --help \n\n")
}

sub open {
my $path = shift;
if ($path !~ m#^/#) {
warnings::warn("changing path to /var/abc")
if warnings::enabled();
$path = "/var/abc/$path";
}
}

sub BetaUpdate {
my $myCurDir = shift;
print "\nNotice: \nYou are going to use your own beta score values\nLooking for the betaScore file at $myCurDir\n"; 
}

sub BetaExists {
print print "Beta score file was found";
}

sub BetaMissed {
my $myCurDir = shift;
print "\nPlease place your own beta score file [ named: betaScore ] in $myCurDir folder.\nFor betaScore file format please check the manual.\n";
exit(0);
}

sub viz {
print "Generating visualization files\t\n";
}

sub ClassifyNoNCBI {
my $myCurDir = shift;
print "\nNotice: \nYou are not using the NCBI taxonomy database to define phylogenetic nodes.\nLooking for the classification.eba file in $myCurDir\n";
}

sub ClassifyCustom {
print "Custom phylogeny file was found";
}

sub yesSPS {
print "File species.sps exists\n"; 
}

sub ClassifyProvide {
my $myCurDir = shift;
print "\nPlease place your classification file [ named: classification.eba ] in $myCurDir folder.\nFor classification.eba file format please check the manual.\n"; 
exit(0);
}

sub EBAWelcome {
my $VERSION = shift;
print "\nWelcome to help section of Evolutionary Breakpoint Analyser (EBA) tool, version $VERSION\n"; 
EBALib::Messages::printUsage(); ## display documentation
}

sub EBAWho {
my $VERSION = shift;
print "\nEvolutionary Breakpoint Analyser, version $VERSION\n"; 
printDeveloper($VERSION);
}

sub EBAMax {
my $VERSION = shift;
print "\nWARNING: \nYou are analyzing more than ten species\nPlease contact the developers for assistance\n";
printDeveloper($VERSION); 
}

sub WorkOnClass {
print "Working on EBR classifications\n"; 
}

sub ParseClass {
print "Parsing the classification file\n";
}

sub ExludeClass {
print "Excluding classification groups defined with only one species\n";
}

sub CalBS {
print "\nCalculating beta scores.. .\t\n";
}

sub ResRange {
my $prime = shift;
print "Your primary resolution $prime provided in the command line is not present. Please fix the error\n"; 
exit(1);
}

sub TimeTaken {
my $startTime = shift;
printf "\nTotal time taken is %.2f minutes\n", (time - $startTime) / 60;
}

sub MergeMSG {
my $curDir=shift;
my $dir=shift;
print "\nMerging resolutions and recalculating probabilities and classifications\n The merged dataset and results are available from $curDir/$dir/Merge\n"; 
}

sub cleanIn {
print "\nCleaning intermediate files\n";
}

sub deleteIn {
print "\nDeleting existing folders\t\n";
}

sub createMSG {
print "Creating the folders required\t\n";
}


sub noop {
my $InFile=shift;
print "open $InFile:";
}

sub noClose {
print "could not close file: $!\n";
}

sub fail {
print "Copy failed:$!";
}

sub failOp {
my $InFile=shift;
print "Could not open $InFile for reading\n";
}

sub failCl {
my $OntFile=shift;
print "Could not close $OntFile\n";
}

sub defEBRs {
print "Defining EBRs\n";
}

sub lookEBR {
my $InFile=shift;
print "Looking for breakpoints in $InFile ...\t\n";
}

# Check the parameters provided by users as commands
sub checkParameters {
my ($options_ref) = @_;
my %options=%$options_ref;
print "Checking parameter\n";
print %options;
# make sure required arguments were passed
my @required = qw( n d r );
}

sub multiHit{
my $InFile=shift;
print "Your $InFile has both scaffold and chromosome entries for the target genome. Please use only one assembly type, either Scaffold or Chromosome\n"; 
exit;
}


sub genomeMSG {
my $InFile=shift;
print "There are wrong chromosome or scaffold values or no entries in column 9 of $InFile\nThere should be chromosome or scaffold ids present in your input file\n"; 
exit(1);
}

sub fileFormat {
print "Warning: file format is not correct.. You might have an extra tab or newline entry at the end of the line\n";
}

sub confused {
my ($InFile, $x, $y, $z)=@_;
print "In $InFile\n HSBs have overlapping coordinates, $x: start $y and end $z.\n To define EBRs, HSBs must be non-overlapping. Please fix your HSB files\n\n";
exit(1); ## Terminate the program 
}

sub build {
my $new=shift;
print "$new EBR table is built...\t\n";
}

sub forgot {
my $chr=shift;
print "WARNING: No chromosome $chr found in your chr_size.txt file\n"; exit(0);
}

sub callProb {
print "Calculating intermediate classification probabilities\t\n";
}

sub callProb2 {
print "Calculating intermediate classification probabilities\nWriting intermediate \t\n";
}

sub noChrOp {
print "Can't open chr_size.txt file\n";
}

sub noChrcl {
print "Could not close file: $!\n";
}

sub betaHai{
print "\nFile beta.tmp exists\n";
}

sub betaDel {
print "-- -> File beta.tmp was deleted\n";
}

sub betaScoreHai {
print "\nFile betaScore exists\n";
}

sub noLineage {
my $fileIn = shift;
print "Please revome lineage name from $fileIn line number $. This will be included automatically\n"; exit(1);
} #did not understand this

sub noGroup {
my $fileIn = shift;
print "Classification group name is missing in $fileIn file at the line number $.\n"; exit(1);
}

sub classError {
my ($species, $tmp1, $file1 , $lineN)= @_;
print die "Species name $species might have a typo or missing values in $tmp1 classification group, please check $file1 file at the line number $lineN\n"; exit(1);
}

sub forgotLineage {
print "Warning: 'lineage=' is missed from the classification.eba file. Please check the manual for a right format\n";
}

sub typoFile {
my ($refName, $fileIn)= @_;
print "The reference genome name $refName, provided in the command line, is missing in $fileIn file. Please check $refName and $fileIn for any typos.\n"; exit(1);
}

sub betaScoreDel {
print "-- - File betaScore deleted\n";
}

sub poissonMSG {
print "\nCalculating final EBR scores and performing classifications\t\n";
}

sub urRes {
my $res= shift;
print "Your resolution is : $res\n";
}

sub noDel {
my $val= shift;
print "Could not remove $val\n";
}

sub numericFold {
my $resolutionName = shift;
print "Folder $resolutionName has not a numeric name and will not be included in the analysis\nHSB folders must have numeric names\n\n"; 
exit(1);
}

sub underscore {
my $file = shift;
print "\nIt seems you HSB file $file name is in a wrong format. Please refer to the manual for more details\t\n"; 
}  

sub duplicateMSG {
my ($dup_ref, $path)=@_;
my @dup=@$dup_ref;
print "\nWe found the duplicate target species name for @dup in $path !! You should not have matching first name of the HSBs files !! \n Please name the file with their species name for better understanding !!\n"; 
exit(1); #can not understand
}

sub noBeta {
print "Beta scores can not be calculated becuase your dataset has less than 3 resolutions. Please build your own betaScore file or provide more resolutions for the analysis\n";  
exit 0;
}

sub zero {
print "Warning: You have 0 missed and 0 real data\n";
} # can not understand what this means

sub tableMSG {
print "Generating the EBR table\t\n";
}

sub spsNum {
my $path=shift;
print "\n The file in $path does not list all expected species\t OR\n The file might not be in the expected format. Please check the manual for more details\t\n"; 
exit(1);
}

sub wrongNum {
my ($InFile, $line)=@_;
print "\nThere are some unexpected values in $InFile file at the line:\n$InFile\n $line\n"; exit(1);
}

sub chrMiss {
my ($InFile, $line)=@_;
print "\nChromosome id is missing in $InFile file at the line:\n$line\n"; exit(1);
}

sub spsMiss {
my ($InFile, $line)=@_;
print "\n Reference and/or target genome name(s) are missing in $InFile file at the line:\n$line\n"; exit(1);
}

sub genomeStatus {
my ($InFile, $line)=@_;
print "\n Target genome assembly status is in a wrong format in $InFile file. It should be either Scaffolds or Chromosomes:\n$line \n"; exit(1);
}

sub noTarSame {
my $InFile  = shift;
print "\n Your HSB file $InFile has more than one unique entry for the target genome name\n"; exit(1);
}

sub noRefSame {
my $InFile  = shift;
print "\nYour HSB file $InFile has more than one unique entry for the reference genome name\n"; exit(1);
}

sub noGenomeSame {
my $InFile  = shift;
print "\nAll the genome information is not same in $InFile file\n It should be either Scaffolds or Chromosomes\n"; exit(1);
}  #can not understand

sub dataVerify {
print "\nVerifying your datasets...\t\n";
}

sub reuseMSG {
print "Calculating the final EBR classification with reuse EBRs...\t\n";
}

sub NA {
my $n=shift;
print "\nWARNING: Line number $n contains some NA or blank spaces\n";
}

sub refId {
my ($inputSpsId, $species)=@_;
print "\nYour reference species taxonomy id and name are: $inputSpsId\t$species\n";
}

sub classProb {
print "Calculating classification probabilities\t\n";
}


sub finalEBR {
print "Generating the final EBR table\t\n";
}

sub classMSG {
print "\n\tThe the taxonomy ids and names of the species in your dataset are:\t\n";
}


sub multiHit {
print "\tWARNING:\n\tMore than one hit was found for the genome name.\nPlease choose any one of them by Number[#]\n\t<#>\t<Ids>\t<Name>\n";
}

sub noHits {
$sp=shift;
print "\t----- FATAL Error -----\n\tNo hits were found for $sp name in our dataset\n";
print "Please provide another name [type a new name in scientific format]:";
}

sub opted {
my ($opted,$k,$sp)=@_;
print "You have opted $opted=$k\t$sp\n";
}

sub printDeveloper {
my $VERSION = shift;
print << "End_Print_Developer";

-------------------------------------------------------------------------------
Script: 
  Evolutionary Breakpoint Framework : Evolutionary Breakpoint Analyser

Aim of the tool:
  Automated definition and classification of evolutionary breakpoints regions from a large number of aligned genomes using species phylogenetic relationships.

Developers: 
  Jitendra Narayan <jin1_at_aber.ac.uk>,  
  Denis Larkin <dmlarkin_at_gmail.com> 

Address:
  Institute of Biological, Environmental & Rural Sciences
  Aberystwyth University 
  Penglais Campus
  Aberystwyth 
  Ceredigion 
  SY23 3DA 

  Tel: +44 (0)1970 628505
  Fax: +44 (0)1970 622350

Version:
  Last updated on 31-01-13 Version: $VERSION

Term:
  This EBA package $VERSION is free software; you can redistribute it and or modify
  it under the terms of the GPL (either version 1, or at your option,
  any later version) or the Artistic License 2.0. 

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
  THE SOFTWARE.
---------------------------------------------------------------------------------

End_Print_Developer
}




# This function prints the script usage mode
sub printUsage {
my ($message) = @_;
if (defined $message) {
    print STDERR "\nERROR: $message\n";
  }

print << "End_Print_Usage";

Usage:

perl EBA.pl -n <number> -d <dir> -r <refname> -p <prime> -t <number> [-o <outdir>]

Mandatory parameters:

  -n <number> 	Provide the number of species that you are going to work with.

  -d <dir> 		The name of directory which contains the list of HSBs. 

  -r <refname> 	Name of your reference species [ Must be scirentific name ].

  -p <prime> 	Provide the primary resolution name [ the resolution name should be numeric ].

  -t <number> 	The threshold value for reuse breakpoint filtration. 

-------------------------------------------------------------------------------
Optional parameters:

  -o <outdir> 	Output directory name (default: EBA_OUT).

  -b  If you would like to calculate the beta then. 

  -c  Want classification of all the studied species.

  -x  If you would like to exclude the classification group with only one species.

  -k  Interested in intermediate files.

  -v  Want to validate your input files.

  -i  <number> Increment the breakpoint sizes by <number> on both side.

  -m To not merge the resolutions

  -l Not to print the lineage informations

---------------------------------------------------------------------------------
File formats:

  + HSBs file:

  The HSBs file has homologous synteny blocks informations, which were identified 
  between two genomes reference(RE) and target(TA) using certain HSB identification softwares ( eq. SyntenyTracker etc.). 
  
  It has the following columns:

  + Column1     - Reference (RE) name
  + Column2     - Reference (RE) chromosome number
  + Column3 	- Start coordinates of RE HSBs
  + Column4   	- End coordinates of the RE HSBs 
  + Column5     - Chromosome name of target(TA) species. 
  + Column6 	- Start coordinates of TA HSBs 
  + Column7		- End coordinates of TA HSBs 
  + Column8		- Orientation of the block 
  + Column9    	- Target (TA) species name
  + Column10	- Target (TA) chromosome [ Scaffolds or Chromosome]

  The breakpoint should must satisfy the following:
     1 : Valid homologous synteny blocks
     2 : Start coordinates of HSBs should be smaller than the end coordiantes
  Warning: The file must NOT have header with the column names. 
  Warinig: The HSBs should not be overlapping.

-------------------------------------------------------------------------------

  The script will generate a FINAL file for each resolutions and write them
  withing each resolution directory. The FINAL file has the following name pattern:
  + Species Names   : Column one to number of species[studied species] columns are the 
		       name of the species and its coordinates with breakpoints and internal score.
  + Chromosome 	    : Chromosome number of the breakpoints.
  + Breakpoints     : Real breakpoints woth both the coodinates seperated by "start--end";
  + Score	     	: Breakpoints score of becoming it as a reference specific. 
  + Narrow Breaks   : Narrow breakpoint sepeated by "stat<->end";
  + Species	     	: Breakpoints and gaps were present in total number of species at this breakpoint.
 
  where: [N] is the number of the breakpoint.

End_Print_Usage

if (defined $message) {
    print "\nERROR: $message\n\n";
  }

exit(1);
}
1;
