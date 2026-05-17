# $Id$ All messages
# Perl module for EBA EBALib::DrawMessages;
# Author: Jitendra Narayan <jnarayan81@gmail.com>
# Copyright (c) 2015 by Jitendra. All rights reserved.
# You may distribute this module under the same terms as Perl itself

##-------------------------------------------------------------------------##
## POD documentation - main docs before the code
##-------------------------------------------------------------------------##

=head1 NAME

EBALib::DrawMessages  - DESCRIPTION of Object

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

sub viz {
print "Generating visualization files\t\n";
}

sub ClassifyCustom {
print "Custom phylogeny file was found";
}

sub yesSPS {
print "File species.sps exists\n"; 
}

sub ParseClass {
print "Parsing the classification file\n";
}

sub ExludeClass {
print "Excluding classification groups defined with only one species\n";
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

sub opted {
my ($opted,$k,$sp)=@_;
print "You have opted $opted=$k\t$sp\n";
}

1;
