# $Id$ Draw Breakpoint Graph
# Perl module for EBA EBALib::Draw::drawBreakpointGraph;
# Author: Jitendra Narayan <jnarayan81@gmail.com>
# Copyright (c) 2015 by Jitendra. All rights reserved.
# You may distribute this module under the same terms as Perl itself

##-------------------------------------------------------------------------##
## POD documentation - main docs before the code
##-------------------------------------------------------------------------##

=head1 NAME

EBALib::Draw::drawBreakpointGraph  - DESCRIPTION of Object

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


package EBALib::Draw::drawBreakpointGraph;

use strict;
use warnings;
#use Term::ANSIColor;
use Cwd;

use Exporter;

our @EXPORT_OK = "breakpointGraph";

sub breakpointGraph {
my ($fileName, $spsArray_ref, $allResolutions_ref, $num)=@_;
my @spsArr=@$spsArray_ref;
my @allRes=@$allResolutions_ref;
my $file2="No_unclassified_EBRs_per_resolution.data";

#print $fileName;
#print @spsArr;
#print @allRes;

my (@allBreaks, $countBreakNumber, @finalData, @all, @catColor);

my @allResolutions = sort { $a <=> $b } @allRes;
my @spsArray = sort { $a cmp $b } @spsArr;
	push (@finalData, \@allResolutions);
	foreach my $speciesName(@spsArray) {    # print "$speciesName\n";
		foreach my $res (@allResolutions) {
			open INFILE,  "$fileName" or die "$0: open $fileName : $!";
				while (<INFILE>) {
					chomp;    
					my $line= lc($_);  #print "$line\n";
					my @tmp = split /\t/, $line;
					s{^\s+|\s+$}{}g foreach @tmp;
					#next if $tmp[3] ne $speciesName;
						if (($tmp[0] == $res) and ($tmp[3] eq $speciesName) and ($tmp[6] eq "break")) {  ## b is small letter
					
							$countBreakNumber++;

						}
				
				}
		
		push (@allBreaks, $countBreakNumber);
		push (@all, $countBreakNumber);
		$countBreakNumber=0;
		}
	my ($rand, $rand2);
	LABEL: { ($rand, $rand2) = random_colors(); redo LABEL if ($rand2 eq "white"); }
	push(@finalData, [@allBreaks]);
	push (@catColor, $rand); #print @catColor;
	undef @allBreaks; 	
	
	
	}
my $maximum=max(@all);
print2d (\@finalData, $file2);
#-generateGraph (\@finalData,\@catColor, \@spsArray, $maximum, $num);
undef @all; undef @finalData; undef @catColor;

}

sub generateGraph {

my ($finalData_ref, $catColor_ref,$allSps_ref, $max, $num)=@_;
my @finalData=@$finalData_ref;
my @catColor=@$catColor_ref;
my @allSps=@$allSps_ref;
my $y_axis_max = $max;

use strict;
use warnings;
#-use GD;
#-use GD::Text;
#-use GD::Graph::bars;


#create graph object for canvas 800 by 600 pixels
my $firstN=(1000+($num*50)); my $secondN=(800+($num*2));

q^
my $my_graph = new GD::Graph::bars($firstN,$secondN, 1);
#set graph options required


$my_graph->set(
	'title' => 'Breakpoints in all resolutions', #graph title
	'y_label' => 'No. EBRs', #y-axis label
	'x_label' => 'HSB definition resolutions',
	'y_max_value' => $y_axis_max, #the max value of the y-axis
	'y_min_value' => 0, #the min value of y-axis, note set below 0 if negative values are required
	'y_tick_number' => 20, #y-axis scale increment
	'y_label_skip' => 2, #label every other y-axis marker
	'box_axis' => 0,
	'line_width' => 2,
	'x_label_position' => .5,
	'y_label_position' => .5, 
	'shadow_depth' => 1,
	'bargroup_spacing' => 20,
	'accent_treshold' => 200, 
	't_margin'=> 10,
	'b_margin'=> 10, 
	'l_margin'=> 10, 
	'r_margin'=> 10,
	#do not draw border around graph
	#width of lines
	'legend_spacing' => 5, #spacing between legend elements
	'legend_placement' =>'RC', #put legend to the centre right of chart
	'dclrs' => \@catColor, #reference to array of category colours for each line
	'transparent' => 0
	) || die "\nFailed to create line graph: $my_graph->error()";
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
#write graph to a file
my $line_file = "No_unclassified_EBRs_per_resolution.gif";

open(IMG, ">$line_file") || die ("\nFailed to save graph to file: $line_file. $!");;
print IMG $plot->gif();
close (IMG);

use Cwd;
my $currentdir= getcwd;
print "Created a line graph for all EBRs in $currentdir/$line_file\n";

q^ if 0;
}



##--------------------------------------------------------------
sub random_colors {
    my ($r, $g, $b) = map { int rand 256 } 1 .. 3;

    my $lum = ($r * 0.3) + ($g * 0.59) + ($b * 0.11);

    my $bg = sprintf("#%02x%02x%02x", $r, $g, $b);
    my $fg = $lum < 128 ? "white" : "black";

    return ($bg, $fg);
}

##-------------------------------------------------------------
sub max {
    my ($max, @vars) = @_;
    for (@vars) {
        $max = $_ if $_ > $max;
    }
    return $max;
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


1;

1;
