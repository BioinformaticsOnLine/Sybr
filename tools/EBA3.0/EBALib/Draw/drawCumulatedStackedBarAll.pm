# $Id$ Draw Cumulative Stack Bar to All
# Perl module for EBA EBALib::Draw::drawCumulatedStackedBarAll;
# Author: Jitendra Narayan <jnarayan81@gmail.com>
# Copyright (c) 2015 by Jitendra. All rights reserved.
# You may distribute this module under the same terms as Perl itself

##-------------------------------------------------------------------------##
## POD documentation - main docs before the code
##-------------------------------------------------------------------------##

=head1 NAME

EBALib::Draw::drawCumulatedStackedBarAll  - DESCRIPTION of Object

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


package EBALib::Draw::drawCumulatedStackedBarAll;

use strict;
use warnings;
#use Term::ANSIColor;

use Exporter;

our @EXPORT_OK = "drawStackBarMergedAll";

sub drawStackBarMergedAll {

use strict;
use warnings;
#-use GD;
#-use GD::Text;
#-use GD::Graph::hbars;

my ($dir, $num)=@_;
my ($InPutFile, $OutPutFile, $OutPutFile2);

$InPutFile="allCumData.tmp";
$OutPutFile="Unique_reuse_EBRs_all_resolutions.gif";
$OutPutFile2="Unique_reuse_EBRs_all_resolutions.data";

##-------------------------------------------------------------

use strict;
use File::Find;

my %myHashNames;my %newHash;  my @allHashValues=();
my @allFinalNames; my @allFinalNames2; my @allFinalValues1; my @allFinalValues2;
find(\&wanted, "$dir");

sub wanted {
#-l && !-e && print "bogus link: $File::Find::name\n";
	if (-f and $File::Find::name =~/Number_unique_reuse_EBRs/) {
	
		#readFile("StackBarGraphFinal.data");
		#print "$File::Find::name\n";
		my @allValues; my $file;
		my @fname = split(/\//, $File::Find::name);
		#print $fname[-1];	
		if ($fname[-1] =~/\.data$/){ 
			$file= $fname[-1]; 
			open(FILE, "$file") || (warn "Can't open \n");
			$/ = "\n";
      			while (<FILE>) { 
				chomp;
				#my @tmp=split(/\t/, $_); 
				push (@allValues, $_);	
				#print "$_\n";
			}
			close FILE;
 		no warnings;
		my @mytmp1=split(/\t/, $allValues[0]);
		my @mytmp2=split(/\t/, $allValues[1]); 
		my @mytmp3=split(/\t/, $allValues[2]);
	
		foreach (my $bb=0; $bb<=$#mytmp1; $bb++) { $myHashNames{$mytmp1[$bb]}="$mytmp2[$bb]:$mytmp3[$bb]"; }
		undef @allValues;
		#print "========\n";
		push (@allHashValues, \%myHashNames);
		#%newHash=%myHashNames; 
		#undef %myHashNames;
    		} 
	}	
}


for(my $xx=0; $xx<=$#allHashValues; $xx++){
my $allHashValues_ref=$allHashValues[$xx];
my %myHashNames= %{$allHashValues_ref}; 
foreach my $key (sort keys %myHashNames) {
		my $val=isInList($key, @allFinalNames);
		next if $val == 1;
		#print "$key\t$myHashNames{$key}\n";
		#%newHash=(%myHashNames,%myHashNames) ;
     		my @tmp=split(/\:/, $myHashNames{$key}); 
     		push @allFinalNames, $key; 
		$key=~ s/_/ /ig;
		push @allFinalNames2, $key; push @allFinalValues1, $tmp[0]; push @allFinalValues2, $tmp[1];
	#print "---------------------------------\n";
	}

#print @allHashValues;
undef %myHashNames;

#print "---------------------------------\n";

}


##--------------------------------------------------------- drawing start here

my @finalData=(\@allFinalNames2, \@allFinalValues1, \@allFinalValues2);

my @allSps=('Unique', 'Reuse'); ## In this script the allSps is actually breaks and pseudobreaks

print2d(\@finalData, $OutPutFile2);

#create graph object for canvas 800 by 600 pixels
my $firstN=(1000+($num*50)); my $secondN=(1000+($num*50));

q^
my $my_graph = new GD::Graph::hbars($firstN,$secondN,1);

#my $my_graph = new GD::Graph::hbars(2000,2000,1);
#set graph options required
#$my_graph->set_legend_font('verdana', 12);
#$my_graph->set_title_font('arial', 20);
$my_graph->set( 
	'title' => 'No.  unique and reuse EBRs per species and phylogenetic group (all resolutions)', #graph title
	'y_label' => 'No.  EBRs', #y-axis label
	#'x_label' => 'Classification_Group',
	'y_min_value' => 0, #the min value of y-axis, note set below 0 if negative values are required
	'y_tick_number' => 20, #y-axis scale increment
	'y_label_skip' => 1, #label every other y-axis marker
	'box_axis' => 1,
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
	'values_space' => 4, 
	'long_ticks' => 1, ## to add lines in background
	#'borderclrs'=> 'black',
	#do not draw border around graph
	#width of lines
	'legend_spacing' => 5, #spacing between legend elements
	#'legend_placement' =>'RC', #put legend to the centre right of chart
	#'dclrs' => \@catColor, #reference to array of category colours for each line
	'dclrs' => [ qw( green lred yellow) ], ## Fix the color ... ## if uncertain added then all three color will work else only 2.
	'cumulate' => 'true', 	
	'transparent' => 0
	
)|| die "\nFailed to create cumulative bar chart: $my_graph->error()"; 

#set legend
$my_graph->set_legend(@allSps); ### This time only Unique and Reuse
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
print "\nCreated a comparative cumulated stacked bar graph for all EBRs at all resolutions in $OutPutFile\n";

q^ if 0;
}

## Subroutines here --------------------------------------


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


1;

