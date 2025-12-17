# $Id$ Concatenate All 
# Perl module for EBA EBALib::ConCatAll;
# Author: Jitendra Narayan <jnarayan81@gmail.com>
# Copyright (c) 2015 by Jitendra. All rights reserved.
# You may distribute this module under the same terms as Perl itself

##-------------------------------------------------------------------------##
## POD documentation - main docs before the code
##-------------------------------------------------------------------------##

=head1 NAME

EBALib::ConCatAll  - DESCRIPTION of Object

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

package EBALib::ConCatAll;
use strict;
use warnings;
#use Term::ANSIColor;

use Exporter;

our @EXPORT_OK = "concatAll";

sub concatAll {
my $path=shift;
open OUTFILE, ">", "$path/EBA_OutFiles/all_all.eba00" or die die $!; 
my $dir="$path/EBA_OutFiles";

opendir(DIR, $dir) or die $!;
   	while (my $file = readdir(DIR)) {
        # We only want files
        next unless (-f "$dir/$file");
        # Use a regular expression to find files ending in .txt
        next unless ($file =~ m/\.eba5$/); 
	my @new_file=split(/_/, $file);
		if ($new_file[0] ne "") { 
			open(FILE, "<$path/EBA_OutFiles/$file") || (warn die EBALib::Messages::failOp($file));
      				while (<FILE>) {
					print OUTFILE $_;
      				}
			close(FILE) or die EBALib::Messages::failCl("FILE");				}				
		  }
closedir(DIR) or die EBALib::Messages::failCl("DIR");
close OUTFILE or die EBALib::Messages::failCl("OUTFILE");
}


1;
