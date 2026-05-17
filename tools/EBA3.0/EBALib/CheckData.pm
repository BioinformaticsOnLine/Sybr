# $Id$ Check Data
# Perl module for EBA EBALib::CheckData;
# Author: Jitendra Narayan <jnarayan81@gmail.com>
# Copyright (c) 2015 by Jitendra. All rights reserved.
# You may distribute this module under the same terms as Perl itself

##-------------------------------------------------------------------------##
## POD documentation - main docs before the code
##-------------------------------------------------------------------------##

=head1 NAME

EBALib::CheckData  - DESCRIPTION of Object

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

package EBALib::CheckData;
use strict;
use warnings;
#use Term::ANSIColor;

use Exporter;

our @EXPORT_OK = "verifyData";

sub verifyData {
EBALib::Messages::dataVerify();

my ($dir, $number)=@_;
my $parent = "./$dir";  #The main folder which contains all the resolution folders
my ($par_dir, $sub_dir);
my @allResolutions;

opendir($par_dir, $parent) or EBALib::Messages::failOp($parent);
while (my $sub_folders = readdir($par_dir)) {
    next if ($sub_folders =~ /^\.\.?$/);  # skip . and ..
    my $path = $parent . '/' . $sub_folders;
    next unless (-d $path);   # skip anything that isn't a directory
	# my @name = split(/\//, $path);
	# print "Verifying $path dataset\n\n";
	loopOver($path, $number);
	my @res = split(/\//, $path);
	my $resolutionName=$res[-1];  ## Name of the resolutions
	if (isInteger($resolutionName)) { }	
	else { EBALib::Messages::numericFold($resolutionName);
	}
} 

#print "\tWell your dataset looks good !!! \n";	
undef @allResolutions;
closedir($par_dir);

} ## loop ends here

sub loopOver {
my ($path, $number)=@_; 
my @speciesNames; my $sub_dir; my @fileNames; my $count=0;
	opendir($sub_dir, $path);
    	while (my $file = readdir($sub_dir)) {
        	next unless $file =~ /\.txt?$/i;
        	my $full_path = $path . '/' . $file;
		my $underscore="_";
			if (index($file, $underscore) == -1) { 
				EBALib::Messages::underscore($file);
				}
		my @first_name = split(/_/, $file); 
		push (@fileNames, $first_name[0]);
			if ($first_name[0] ne "")  {
				my $InFileName="$path/$file";
        			fileValidation ($path, $InFileName);
			}
	$count++; 	
	}
	my %seen = (); my @dup = map { 1==$seen{$_}++ ? $_ : () } @fileNames; 
	if(@dup) { EBALib::Messages::duplicateMSG(\@dup, $path);}
	undef @fileNames; undef %seen;
	#print "$count\n"; ## Checking for number of files in your folders..
	if ($count != $number) { 
		EBALib::Messages::spsNum($path);
		}
	closedir($sub_dir);
}

sub fileValidation {
my ($path, $InFile)=@_;
open INFILE,  $InFile or die EBALib::Messages::failOp($InFile);
$/ = "\n";
my $decision; my @allreference; my @alltarget; my @ref; my @tar; my @genomeInfo;my @info; my @allChr; my @chr;
while (<INFILE>) {
	my $line=lc($_);
	chomp $line;
	my @tmp=split /\t/, $line;
	for (@tmp)  { s/^\s+//; s/\s+$//; }#replace one or more spaces 
		if (isInteger($tmp[2]) and isInteger($tmp[3]) and isInteger($tmp[5]) and isInteger($tmp[6])) { }
		else { EBALib::Messages::wrongNum($InFile, $line);}
		
		if(($tmp[1]=~/^\s*$/) or ($tmp[4]=~/^\s*$/)){ EBALib::Messages::chrMiss($InFile, $line);}
		#if (isSign($tmp[7])) { } else { print " The block orientation is not provided in this line: \n$line\n \n"; exit(1); } 
		if (($tmp[0]=~/^\s*$/) or ($tmp[8]=~/^\s*$/)){ EBALib::Messages::spsMiss($InFile, $line);}
		if (($tmp[9] ne "chromosomes") and ($tmp[9] ne "scaffolds")) { EBALib::Messages::genomeStatus($InFile, $line);}	
		
	push(@allreference,$tmp[0]);
	push(@alltarget,$tmp[8]);
	push(@genomeInfo,$tmp[9]);
	push(@allChr, $tmp[1]);

	}
@ref = EBALib::CommonSubs::getUniq(@allreference);
@tar = EBALib::CommonSubs::getUniq(@alltarget);
@info = EBALib::CommonSubs::getUniq(@genomeInfo);
@chr=  EBALib::CommonSubs::getUniq(@allChr);

my @chrAll;
open(CHRFILE1, "chr_size.txt") || warn EBALib::Messages::failOp("chr_size");
while (<CHRFILE1>) { chomp $_; next if $_ =~ /^$/; next if $_ =~ /^\s*#/; my @chrtmp = split /\t/, lc($_); push @chrAll,$chrtmp[0];}       ## We need to improve it !!!!!!!!!1
close CHRFILE1 or die EBALib::Messages::failCl("file");
for (@chrAll) { s/^\s+//; s/\s+$//;} # replace spaces


if($#tar > 0) { EBALib::Messages::noTarSame($InFile);}
elsif ($#ref > 0) { EBALib::Messages::noRefSame($InFile);}
elsif ($#info > 0) { EBALib::Messages::noGenomeSame($InFile);}
#elsif ($#chr != $#chrAll) { print "It looks you have wrong chromosome size file where the number of chromosome is different from that of your reference species\n"; exit(1);} 

undef @allreference;
undef @alltarget;
undef @genomeInfo;
undef @allChr;
close INFILE or die EBALib::Messages::failCl($InFile);
}

## ---------------------------------------------------------------------------- 
#subroutines here

sub isInteger { defined $_[0] && $_[0] =~ /^\d+$/; }
sub isSign { defined $_[0] && $_[0] =~ /^[+-]$/; } 



1;
