# $Id$ Draw Beta Graph
# Perl module for EBA EBALib::Draw::drawBetaGraph;
# Author: Jitendra Narayan <jnarayan81@gmail.com>
# Copyright (c) 2015 by Jitendra. All rights reserved.
# You may distribute this module under the same terms as Perl itself

##-------------------------------------------------------------------------##
## POD documentation - main docs before the code
##-------------------------------------------------------------------------##

=head1 NAME

EBALib::Draw::drawBetaGraph  - DESCRIPTION of Object

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

package EBALib::Draw::drawBetaGraph;

use strict;
use warnings;
#-use GD;
#-use GD::Text;
#-use GD::Graph::lines;
#use Term::ANSIColor;
#-use GD::Graph::mixed;

use Exporter;

our @EXPORT_OK = "drawBetaGraph";

sub drawBetaGraph { # for line graph

my ($finalData_ref, $catColor_ref,$allSps_ref, $max)=readData('betaScore');
my @finalData=@$finalData_ref;
my @catColor=@$catColor_ref;
my @allSps=@$allSps_ref;
my $y_axis_max = $max;
my $file2='betaScore.data';
#create graph object for canvas 800 by 600 pixels

q^
my $my_graph = new GD::Graph::mixed(1000,600,1);
#set graph options required

$my_graph->set(
'title' => 'Beta Scores in all resolutions', #graph title
'y_label' => 'Beta Scores', #y-axis label
'x_label' => 'Resolutions',
'y_max_value' => $y_axis_max, #the max value of the y-axis
'y_min_value' => 0, #the min value of y-axis, note set below 0 if negative values are required
'y_tick_number' => 20, #y-axis scale increment
'y_label_skip' => 2, #label every other y-axis marker
'box_axis' => 0,
'line_width' => 2,
'x_label_position' => .5,
'y_label_position' => .5, 
#do not draw border around graph
#width of lines
'legend_spacing' => 5, #spacing between legend elements
't_margin'=> 20,
'b_margin'=> 20, 
'l_margin'=> 20, 
'r_margin'=> 20,
'legend_placement' =>'RC', #put legend to the centre right of chart
'dclrs' => \@catColor, #reference to array of category colours for each line
'transparent' => 0
) || die die "\nFailed to create line graph: $my_graph->error()";

#set legend
$my_graph->set_legend(@allSps);
#plot graph with table data
$my_graph->set_title_font(['verdana', 'arial', gdMediumBoldFont],14); 
$my_graph->set_x_label_font(['verdana', 'arial', gdMediumBoldFont],12); 
$my_graph->set_y_label_font(['verdana', 'arial', gdMediumBoldFont],12); 
$my_graph->set_x_axis_font(['verdana', 'arial', gdMediumBoldFont],12); 
$my_graph->set_y_axis_font(['verdana', 'arial', gdMediumBoldFont],12); 
$my_graph->set_legend_font(['verdana', 'arial', gdMediumBoldFont],12); 

my $plot = $my_graph->plot(\@finalData);

q^ if 0;
print2d (\@finalData, $file2);
#write graph to a file

q^
my $line_file = "betaScore.gif";

open(IMG, ">$line_file") || die die ("\nFailed to save graph to file: $line_file. $!");;
print IMG $plot->gif();
close (IMG);

use Cwd;
my $currentdir= getcwd;
print "Created a line graph for beta scores in $currentdir/$line_file\n";
q^ if 0;
}

## Subroutines here --------------------------------------

sub readData {
my $file=shift;
my @allResolutions; my @allValues; my @allSpecies;
open(FILE, "$file") || (warn die "Can't open file $_\n");
      while (<FILE>) { 
	chomp;
	my @tmp=split(/\t/, $_); 
	my @resNname=split(/\:/,$tmp[0]);
	push (@allResolutions,$resNname[0]);
	push (@allValues, $_);
	push (@allSpecies, $resNname[1]);
	print "$resNname[0]\t$resNname[1]\t$_\n";	
	}	
close(FILE);

my @uniqRes=uniq(@allResolutions);
my @uniqResolutions = sort { $a <=> $b } @uniqRes;
my @uniqSps=uniq(@allSpecies);
my @uniqSpecies = sort { $a cmp $b } @uniqSps;
my @finalData; my @catColor; my @scores; my @all;

push (@finalData, \@uniqResolutions);
#push (@finalData, \@uniqSpecies);

	foreach my $sps(@uniqSpecies) { 
	#foreach my $res(@uniqResolutions) {
		foreach my $res(@uniqResolutions) {
		#foreach my $sps(@uniqSpecies) {  
			foreach my $val (@allValues) { 
				my @line=split(/\t/, $val); #print "$line[1]\n";
				my @nr=split(/\:/,$line[0]);
			
				if (($nr[0] == $res) and ($nr[1] eq $sps))  { 
					my $newVal=$line[1]*100;
					my $newNum=sprintf("%.2f", $newVal);
					push (@scores, $newNum);
					push (@all, $newNum);
					#print "$res:$line[1]\n";
				}
			}
		}
	my ($rand, $rand2) = random_colors();

	#push (@doneColor, $rand);
	# if ( grep( /^$value$/, @array ) ) { } else {}  ## We need to check all the done color to avoide matching color 

	push(@finalData, [@scores]);
	push (@catColor, $rand); #print @catColor;
	undef @scores; 
	}	
my $maximum=max(@all);
return (\@finalData,\@catColor, \@uniqSpecies, $maximum);
undef @all; undef @finalData; undef @catColor; undef @allResolutions; undef @allValues; undef @allSpecies;
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
