# $Id$ Check Classification
# Perl module for EBA EBALib::CheckClassification;
# Author: Jitendra Narayan <jnarayan81@gmail.com>
# Copyright (c) 2015 by Jitendra. All rights reserved.
# You may distribute this module under the same terms as Perl itself

##-------------------------------------------------------------------------##
## POD documentation - main docs before the code
##-------------------------------------------------------------------------##

=head1 NAME

EBALib::CheckClassification  - DESCRIPTION of Object

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

package EBALib::CheckClassification;
use strict; 
use warnings;
#use Term::ANSIColor;

use Exporter;

our @EXPORT_OK = "verifyClassification";

sub verifyClassification {
    my ($fileIn, $refName) = @_;
    open my $fh, '<', $fileIn or die EBALib::Messages::failOp($fileIn);
	
    my @spstmp;
    open(SPSFILE, "sps.txt") || warn EBALib::Messages::failOp("file");
    while (<SPSFILE>) { chomp; next if $_ =~ /^\s*#/; @spstmp = split /\,/, lc($_);}
    close SPSFILE or die EBALib::Messages::failCl("file");

    my @allClassification; my @speciesNames;
    while (<$fh>) {
	chomp $_; my $line=lc($_);
	next if $line =~ /^\s*$/;
	if (index($line,"#") == 0) { next; } # Lines starting with a hash mark are comments
	$line=EBALib::CommonSubs::trim($line);
	$line=~ s/\r//g;
	next if $line =~ /^$/; ##blank line
	my @tmpLine=split /\=/, $line; 
	if ($line =~ /^lineage/) { if (defined($tmpLine[1])) { EBALib::Messages::noLineage($fileIn)} } else { @speciesNames=split /\,/, $tmpLine[1]; }

	if ($tmpLine[0] =~/^\s*$/) {EBALib::Messages::noGroup($fileIn) } else { push @allClassification, lc($tmpLine[0]);}

	foreach my $species (@speciesNames) { my $res=EBALib::CommonSubs::isInList($species, @spstmp); if ($res == 0) { EBALib::Messages::classError($species, $tmpLine[0], $fileIn , $.)}}
	undef @speciesNames;
    }
my $group=EBALib::CommonSubs::isInList('lineage', @allClassification);
if($group == 0) { EBALib::Messages::forgotLineage();}
my $refCheck=EBALib::CommonSubs::isInList(lc($refName), @allClassification);
if($refCheck == 0) { EBALib::Messages::typoFile($refName, $fileIn);}
undef @allClassification; 
}

1;



