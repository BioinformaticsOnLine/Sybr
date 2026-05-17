# $Id$ All messages
# Perl module for EBA EBALib::CommonSubs;
# Author: Jitendra Narayan <jnarayan81@gmail.com>
# Copyright (c) 2015 by Jitendra. All rights reserved.
# You may distribute this module under the same terms as Perl itself

##-------------------------------------------------------------------------##
## POD documentation - main docs before the code
##-------------------------------------------------------------------------##

=head1 NAME

EBALib::CommonSubs  - DESCRIPTION of Object

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

package EBALib::CommonSubs;
use warnings::register;

# Checks if a provided two coordinates overlaps or not
sub checkCorOverlaps {
my ($x1, $x2, $y1, $y2)=@_;
return $x1 <= $y2 && $y1 <= $x2;
}

sub open {
my $path = shift;
if ($path !~ m#^/#) {
warnings::warn("changing path to /var/abc")
if warnings::enabled();
$path = "/var/abc/$path";
}
}

# Checks if a provided element exists in the provided list
# Usage: isInList <needle element> <haystack list>
# Returns: 0/1
sub isInList {
	my $needle = shift;
	my @haystack = @_;
	foreach my $hay (@haystack) {
		if ( $needle eq $hay ) {
			return 1;
		}
	}
	return 0;
}

## To find the unique
sub uniq {
    my %seen = ();
    my @return = ();
    foreach my $value (@_) {
        unless ($seen{$value}) {
            push @return, $value;
            $seen{$value} = 1;
        }
    }
    return @return;
}

#sub uniq { my %seen; return grep { !$seen{$_}++ } @_; }

# Perl trim function to remove whitespace from the start and end of the string
sub trim($)
{
	my $string = shift;
	$string =~ s/^[\t\s]+//;
	$string =~ s/[\t\s]+$//;
	$string =~ s/[\r\n]+$//; ## remove odd or bad newline ...
	return $string;
}

# Get the unique values in an array
sub getUniq {
my @set = @_;
my %seen;
my @unique = grep { not $seen{$_} ++ } @set;
return @unique;
undef %seen;
}

# Find the index values in an array
sub indexArray($@) {  
 	my $s=shift;
 	$_ eq $s && return @_ while $_=pop;
 	-1;
}

# To sort unique hash 
sub sortUniqueHash {
    my %hash;
    @hash{@_} = ();
    return sort keys %hash;
undef %hash;
}

sub isInteger { defined $_[0] && $_[0] =~ /^\d+$/; }

# To print any hash 
sub printhash {
  	my %hash=%{$_[0]};
	foreach my $key (sort keys %hash) {
     	print "$key : $hash{$key}\n";
	}
}

# To print any 2d array 
sub print2d {
	my @array_2d=@_;
	for(my $i = 0; $i <= $#array_2d; $i++){
	   for(my $j = 0; $j <= $#{$array_2d[0]} ; $j++){
	      print "$array_2d[$i][$j]\t";
	   }
	   print "\n";
	}
}

## Subroutines NOT appied in EBA
sub present_in {
my $val;
my @a = @{ $_[0] };  # remember to use "my"!
my @b = @{ $_[1] };  # remember to use "my"!

# print "@a\n";
my %seen = ();
my @union = grep { $seen{ $_ }++ } @a, @b;

if (@union) {
	if ($#union eq $#a) { $val= 1;} else { $val= 0;}
	#print "$_\n" for @union;
	}
return $val;  undef %seen;
}


# To expand the scientific name like 10E-2
sub expand {
        my $n = shift;
        return $n unless $n =~ /^(.*)e([-+]?)(.*)$/;
        my ($num, $sign, $exp) = ($1, $2, $3);
        my $sig = $sign eq '-' ? "." . ($exp - 1 + length $num) : '';
        return sprintf "%${sig}f", $n;
}

# To print any 2d array 
sub print2dOTHER {
my ($array2d_ref)=@_;
my @array2d=@$array2d_ref;

foreach my $row(@array2d){
   foreach my $val(@$row){
      print "$val\t";
   }
   print "\n";
}

}


# Find min max in an array
sub MinMax {
my @array=@_;
my ($min, $max);
if(scalar(@array) == 1) {$min=$array[0];} elsif (scalar(@array) == 0) { $min=0;}
else {
	for (@array) {
    		$min = $_ if !$min || $_ < $min;
    		$max = $_ if !$max || $_ > $max
	}
     }
return $min;
}


# Perl doesn't have round, so let's implement it
sub round { my($number) = shift; return int($number + .5 * ($number <=> 0)); }


# To multiply all values in an array with each other
sub multiply_all {
my @all_values = @_;
my $score=1;
foreach my $xx (@all_values) { #print "$xx\t------------\n"; 
	$score = $xx * $score;}
	if (scalar (@all_values) == 0) { $score=0;} 
return $score;
}

# Max from an array calculation
sub max {
    my ($max, @vars) = @_;
    for (@vars) {
        $max = $_ if $_ > $max;
    }
    return $max;
}


## Pause the program, activate when press enter
sub pauseProgram($) {
	my $sleep = shift;    #force sleep;
	if ($sleep) { sleep($sleep); }
	my $flag = 0;
	print "[press enter] : ";
	my $input = '';
	while ( !$input ) {
		print "\b";
		$input = <STDIN>;
	}
}

## returns 1 for y
sub yesORno() { 
	my $yflag = 0;
	print "[y/n] : ";
	my $input = '';
	while ( $input !~ /y|n/i ) {
		print "\b";
		$input = <STDIN>;
		chomp $input;
	}
	if ( $input =~ /^y/i ) { $yflag = 1; }
	return $yflag;
}

sub dircopy {
my @dirlist=($_[0]);
my @dircopy=($_[1]);
until (scalar(@dirlist)==0) {
	mkdir "$dircopy[0]";
	opendir my($dh),$dirlist[0];
	my @filelist=grep {!/^\.\.?$/} readdir $dh;
	for my $i (0..scalar(@filelist)-1) {
		if ( -f "$dirlist[0]/$filelist[$i]" ) {
			EBALib::CommonSubs::fcopy("$dirlist[0]/$filelist[$i]","$dircopy[0]/$filelist[$i]");
		}
		if ( -d "$dirlist[0]/$filelist[$i]" ) {
			push @dirlist,"$dirlist[0]/$filelist[$i]";
			push @dircopy,"$dircopy[0]/$filelist[$i]";
		}
	}
	closedir $dh;
	shift @dirlist;shift @dircopy;
}
}

sub fcopy {
my ($i,$data,$cpo,$cpn);
open($cpo,"<",$_[0]) or die $!; binmode($cpo);
open($cpn,">",$_[1]) or die $!; binmode($cpn);
while (($i=sysread $cpo,$data,4096)!=0){print $cpn $data};
close($cpn);close($cpo);
}



sub returnOS($) {
use Config;
my $EBAv = shift;
if (my $distro = $Config{osname}) {
      my $version = $Config{archname};
      print "\nYou are running EBA version $EBAv on your OS $distro, version $version ---\n";
            if ($distro ne "linux") {
		print "This EBA $EBAv is not supported on this $distro OS; Run at your own risk\n";
	    }
} 
else {
      print "I am afraid, can't recognise your OS distribution\n";
}

}


1;
