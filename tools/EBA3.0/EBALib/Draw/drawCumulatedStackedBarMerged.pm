# $Id$ Draw Cumulative Bar Merged
# Perl module for EBA EBALib::Draw::drawCumulatedStackedBarMerged;
# Author: Jitendra Narayan <jnarayan81@gmail.com>
# Copyright (c) 2015 by Jitendra. All rights reserved.
# You may distribute this module under the same terms as Perl itself

##-------------------------------------------------------------------------##
## POD documentation - main docs before the code
##-------------------------------------------------------------------------##

=head1 NAME

EBALib::Draw::drawCumulatedStackedBarMerged - DESCRIPTION of Object

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

package EBALib::Draw::drawCumulatedStackedBarMerged;

use strict;
use warnings;
#use Term::ANSIColor;

use Exporter;

our @EXPORT_OK = "drawStackBarMerged";

sub drawStackBarMerged {

my ($num)=@_;
use strict;
use warnings;
#-use GD;
#-use GD::Text;
#-use GD::Graph::hbars;

my ($InPutFile, $OutPutFile, $OutPutFile2);

$InPutFile="Result_Reuse_Merge.final";
$OutPutFile="Number_unique_reuse_EBRs.gif";
$OutPutFile2="Number_unique_reuse_EBRs.data";

my ($finalData_ref, $catColor_ref, $allSps_ref)=readData($InPutFile, $OutPutFile2, "Merge");
my @finalData=@$finalData_ref;
my @catColor=@$catColor_ref;
my @allSps=@$allSps_ref; ## In this script the allSps is actually breaks and pseudobreaks

#create graph object for canvas 800 by 600 pixels
my $firstN=(800+($num*5)); my $secondN=(800+($num*5));
q^
my $my_graph = new GD::Graph::hbars($firstN,$secondN, 1);
#set graph options required
#$my_graph->set_legend_font('verdana', 12);
#$my_graph->set_title_font('arial', 20);
$my_graph->set( 
	'title' => 'No. unique and reuse EBRs per species and phylogenetic group', #graph title
	'y_label' => 'No. EBRs', #y-axis label
	#'x_label' => 'Classification_Group',
	'y_min_value' => 0, #the min value of y-axis, note set below 0 if negative values are required
	'y_tick_number' => 10, #y-axis scale increment
	'y_label_skip' => 1, #label every other y-axis marker
	#'box_axis' => 1,
	'line_width' => 2,
	'x_label_position' => .5,
	'y_label_position' => .5, 
	'shadow_depth' => 1,
	'bargroup_spacing' => 4,
	'accent_treshold' => 200, 
	't_margin'=> 10,
	'b_margin'=> 10, 
	'l_margin'=> 10, 
	'r_margin'=> 10,
	#'show_values' => 1,
	#'long_ticks' => 0, ## to add lines in background
	'y_long_ticks' => 1,
	'values_space' => 4, 
	#'borderclrs'=> 'black',
	#do not draw border around graph
	#width of lines
	'legend_spacing' => 5, #spacing between legend elements
	#'legend_placement' =>'RC', #put legend to the centre right of chart
	#'dclrs' => \@catColor, #reference to array of category colours for each line
	'dclrs' => [ qw( green lred yellow) ], ## Fix the color ... ## if uncertain added then all trhee color will work else only 2.
	'cumulate' => 'true', 	
	'transparent' => 0
	
)|| die "\nFailed to create cumulative bar chart: $my_graph->error()"; 

#set legend
$my_graph->set_legend(@allSps);
#$my_graph->set_legend_font('verdana', 2); ## not working right now ...
#plot graph with table data

$my_graph->set_title_font(['verdana', 'arial', gdMediumBoldFont],16); 
$my_graph->set_x_label_font(['verdana', 'arial', gdMediumBoldFont],12); 
$my_graph->set_y_label_font(['verdana', 'arial', gdMediumBoldFont],12); 
$my_graph->set_x_axis_font(['verdana', 'arial', gdMediumBoldFont],12); 
$my_graph->set_y_axis_font(['verdana', 'arial', gdMediumBoldFont],12); 
$my_graph->set_legend_font(['verdana', 'arial', gdMediumBoldFont],12);

my $plot = $my_graph->plot(\@finalData);

open(IMG, ">$OutPutFile") || die ("\nFailed to save graph to file: $OutPutFile. $!");
print IMG $plot->gif();
close (IMG);

use Cwd;
my $currentdir= getcwd;
print "\nCreated a cumulated stacked bar graph for merged final EBRs in $currentdir/$OutPutFile\n";

q^ if 0;
}

## Subroutines here --------------------------------------

sub readData {
my ($file, $file2, $resN)= @_;
my @allResolutions; my @allValues; my @allSpecies; my @decision;
open(FILE, "$file") || (warn "Can't open file $_\n");
      while (<FILE>) { 
	chomp;
	next if $. == 1;
	##if header file .. need to next
	my @tmp=split(/\t/, $_); 
	next if $tmp[5] eq "Uncertain"; ## No considering uncertain case
	push (@allValues, $_);
	my @nam=split(/\:/, $tmp[2]); 
	push (@allSpecies, "$nam[0]::$resN");
	push (@decision, $tmp[5]);
	
	}
close(FILE);

my @uniqSps=uniq(@allSpecies);
my @decisionUniq=uniq(@decision);
my @uniqSpecies = sort { $a cmp $b } @uniqSps;
my @finalData; my @finalData2; my @catColor; my @scores; my @all; my @uniqSpecies2;

foreach (@uniqSpecies) { my @mSps=split(/\:\:/, $_); $mSps[0]=~ s/_/ /ig; push (@uniqSpecies2, $mSps[0]);}
push (@finalData2, \@uniqSpecies2);

push (@finalData, \@uniqSpecies);
	foreach my $decision (@decisionUniq){
	foreach my $sps(@uniqSpecies) { 
		my $counter=0; 
		foreach my $val (@allValues) { 
			my @line=split(/\t/, $val);	
			my @name=split(/\:/, $line[2]);
			my @mySps=split(/\:\:/, $sps);
			if (($mySps[0] eq $name[0]) and ($line[5] eq $decision)) { $counter++; }
		}
	if($counter == 0) { $counter=undef;}
	push (@scores, $counter); 
	$counter=0; 
	}
	my ($rand, $rand2);
	LABEL: { ($rand, $rand2) = random_colors(); redo LABEL if ($rand2 eq "white"); }
	push (@catColor, $rand); #print @catColor;
	push(@finalData, [@scores]);
	push(@finalData2, [@scores]);	
	undef @scores;
}
print2d (\@finalData, $file2);
return (\@finalData2,\@catColor, \@decisionUniq); ## provided the breaks and pseudobreaks
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
	no warnings;
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

