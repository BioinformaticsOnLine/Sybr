# $Id$ Draw Final Breakpoint Pie Chart for Merged
# Perl module for EBA EBALib::Draw::drawPieChartBreaksFinalMerged;
# Author: Jitendra Narayan <jnarayan81@gmail.com>
# Copyright (c) 2015 by Jitendra. All rights reserved.
# You may distribute this module under the same terms as Perl itself

##-------------------------------------------------------------------------##
## POD documentation - main docs before the code
##-------------------------------------------------------------------------##

=head1 NAME

EBALib::Draw::drawPieChartBreaksFinalMerged - DESCRIPTION of Object

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


package EBALib::Draw::drawPieChartBreaksFinalMerged;

use strict;
use warnings;
#use Term::ANSIColor;

use Exporter;

our @EXPORT_OK = "drawBreaksPieMerged";

sub drawBreaksPieMerged {
my ($spsNum) = @_;


use strict;
use warnings;
#-use GD;
#-use GD::Text;
#-use GD::Graph::pie; 

my ($InPutFile, $OutPutFile, $OutPutFile2);

$InPutFile="Result_Reuse_Merge.final";
$OutPutFile="EBR_classification_fractions.gif";
$OutPutFile2="EBR_classification_fractions.data";

my ($finalData_ref, $catColor_ref)=readData($InPutFile, $OutPutFile2);
my @finalData=@$finalData_ref;
my @catColor=@$catColor_ref;

#create graph object for canvas 800 by 600 pixels
my $firstN=(800+(2*$spsNum)); my $secondN=(800+(2*$spsNum));

q^
my $my_graph = new GD::Graph::pie($firstN,$secondN, 1);
#set graph options required

$my_graph->set( 
'title'=> "Fraction of unique, reuse and uncertain EBRs at merge resolution (percentage of EBRs)",
'axislabelclr' => 'black', #colour of label segements
'accentclr' => 'black', #colour dividing segements
'start_angle' => 90, 
'3d' => 0, 
#'label' => 'Final Breaks Pie Chart',
# The following should prevent the 7th slice from getting a label 
'suppress_angle' => 5, 
'pie_height' => 36,
't_margin'=> 20,
'b_margin'=> 20, 
'l_margin'=> 20, 
'r_margin'=> 20, 
'transparent' => 0, 
#'dclrs'=>\@catColor,
'dclrs' =>[ qw( green lred cyan) ], ## Fix the color ...
)|| die "\nFailed to create pie chart: $my_graph->error()"; 

$my_graph->set_title_font(['verdana', 'arial', gdMediumBoldFont],16); 
$my_graph->set_label_font(['verdana', 'arial', gdMediumBoldFont],12); 
$my_graph->set_label_font(['verdana', 'arial', gdMediumBoldFont],12);  
$my_graph->set_value_font(['verdana', 'arial', gdMediumBoldFont],12);

my $plot = $my_graph->plot(\@finalData);

open(IMG, ">$OutPutFile") || die ("\nFailed to save graph to file: $OutPutFile. $!");
print IMG $plot->gif();
close (IMG);

use Cwd;
my $currentdir= getcwd;
print "\nCreated a pie chart for merged final EBRs in $currentdir/$OutPutFile\n";
q^ if 0;
}

## Subroutines here --------------------------------------

sub readData {
use Math::Round;
my ($file, $file2, $spsNum)= @_;
my @allResolutions; my @allValues; my @allClass;
open(FILE, "$file") || (warn "Can't open file $_\n");
      while (<FILE>) { 
	chomp;
	if ($_ =~ /^\s*#/) { next; }
	next if $. == 1;
	my @tmp=split(/\t/, $_); 
	push (@allValues, $_);

	push (@allClass, $tmp[5]);	
	}	
close(FILE);

my @uniqClass=uniq(@allClass);
my @uniqClassSorted = sort { $a cmp $b } @uniqClass;
my @uniqClassSortedFinal= grep {$_} @uniqClassSorted; ## Remove the empty index
my @finalData; my @catColor; my @scores; my @all;
my @newUniqSpsName;

#push (@finalData, \@uniqClassSortedFinal);

	foreach my $sps(@uniqClassSortedFinal) { 
		my $counter=0;
		foreach my $val (@allValues) { 
			my @line=split(/\t/, $val);
			if ($sps eq $line[5]) { $counter++; }
		}
	push (@scores, $counter);
	$sps=~ s/_/ /ig;  
	my $ScorePer=round (($counter/scalar(@allValues))*100);
	my $str=join(" ", $sps,"($ScorePer%)");
	push (@newUniqSpsName, $str);
	my ($rand, $rand2);
	LABEL: { ($rand, $rand2) = random_colors(); redo LABEL if ($rand2 eq "white"); }
	push (@catColor, $rand); #print @catColor;
	$counter=0; 
	}
	push(@finalData, \@newUniqSpsName);
	push(@finalData, [@scores]);	
	undef @scores;

print2d (\@finalData, $file2);
return (\@finalData,\@catColor, \@newUniqSpsName);
}

##-------------------------------------------------------------
sub max {
    my ($max, @vars) = @_;
    for (@vars) {
        $max = $_ if $_ > $max;
    }
    return $max;
}

##--------------------------------------------------------------
sub random_colors {
    my ($r, $g, $b) = map { int rand 256 } 1 .. 3;

    my $lum = ($r * 0.3) + ($g * 0.59) + ($b * 0.11);

    my $bg = sprintf("#%02x%02x%02x", $r, $g, $b);
    my $fg = $lum < 128 ? "white" : "black";

    return ($bg, $fg);
}

##----------------------------------------------------------------
## to print any 2d array 
sub print2d {
my ($array2d_ref, $file2)=@_;
my @array2d=@$array2d_ref;

open OUTFILE, ">" , $file2 or die "$0: open $file2: $!";

#print "Print Using ForEach\n";
foreach my $row(@array2d){
   foreach my $val(@$row){
      print OUTFILE "$val\t";
   }
   print OUTFILE "\n";
}
close OUTFILE or die "could not close file: $!\n";
}

##-----------------------------------------------------------------
sub uniq { return keys %{{ map { $_ => 1 } @_ }};}

##------------------------------------------------------------------

1;

