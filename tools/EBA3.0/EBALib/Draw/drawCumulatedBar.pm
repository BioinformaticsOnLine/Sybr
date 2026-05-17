# $Id$ Draw Cumulative Bar 
# Perl module for EBA EBALib::Draw::drawCumulatedBar;
# Author: Jitendra Narayan <jnarayan81@gmail.com>
# Copyright (c) 2015 by Jitendra. All rights reserved.
# You may distribute this module under the same terms as Perl itself

##-------------------------------------------------------------------------##
## POD documentation - main docs before the code
##-------------------------------------------------------------------------##

=head1 NAME

EBALib::Draw::drawCumulatedBar  - DESCRIPTION of Object

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

package EBALib::Draw::drawCumulatedBar;

use strict;
use warnings;
#use Term::ANSIColor;

use Exporter;

our @EXPORT_OK = "drawCumBar";

sub drawCumBar {
my ($path, $resName, $num) = @_;
my $dir="$path/EBA_OutFiles";
print "Making a bar chart for EBRs and pseudobreaks for each genome.. .\t\n";

		 opendir(DIR, $dir) or die $!;
   		 while (my $file = readdir(DIR)) {
        		# We only want files
        		next unless (-f "$dir/$file");
        		# Use a regular expression to find files ending in .txt
        		next unless ($file =~ m/\.eba0$/); 
			my @new_file=split(/_/, $file);
				if ($new_file[0] ne "") { 
					my $InPutFile="$path/EBA_OutFiles/all_brk.eba0";
					my $OutPutFile= "$path/EBA_ImageFiles/Fraction_EBRs_and_gaps_$resName.gif";
					my $OutPutFile2= "$path/EBA_ImageFiles/Fraction_EBRs_and_gaps_$resName.data";
					drawCumBarChart($InPutFile, $OutPutFile, $OutPutFile2, $resName, $num);

				}				
		  }
		closedir(DIR) or die "Can't close DIR: $!";
}

sub drawCumBarChart {

use strict;
use warnings;
#-use GD;
#-use GD::Text;
#-use GD::Graph::hbars;

my ($InPutFile, $OutPutFile, $OutPutFile2, $resName, $num) = @_;

my ($finalData_ref, $catColor_ref, $allSps_ref)=readData($InPutFile, $OutPutFile2);
my @finalData=@$finalData_ref;
my @catColor=@$catColor_ref;
my @allSps=@$allSps_ref; ## In this script the allSps is actually breaks and pseudobreaks

#create graph object for canvas 800 by 600 pixels
my $firstN=(800+($num*5)); my $secondN=(800+($num*5));

q^
my $my_graph = new GD::Graph::hbars($firstN,$secondN, 1);
#set graph options required

$my_graph->set( 
	'title' => "No. unclassified EBRs and in between scaffold gaps ($resName Kbp)", #graph title $resName Kb resolution
	'y_label' => 'No. EBRs and gaps', #y-axis label
	#'x_label' => 'Species_Names',
	'y_min_value' => 0, #the min value of y-axis, note set below 0 if negative values are required
	'y_tick_number' => 10, #y-axis scale increment
	'y_label_skip' => 1, #label every other y-axis marker
	'box_axis' => 1,
	'line_width' => 2,
	'x_label_position' => .5,
	'y_label_position' => .5, 
	'shadow_depth' => 1,
	'bargroup_spacing' => 4,
	'accent_treshold' => 200, 
	't_margin'=> 20,
	'b_margin'=> 20, 
	'l_margin'=> 20, 
	'r_margin'=> 20,
	#'long_ticks' => 1, ## to add lines in background
	#'x_long_ticks'=> 4,
	'y_long_ticks'=> 4,
	#'x_tick_length'=> 4,
	#'y_tick_length'=> 4,
	#'show_values' => 1,
	'values_space' => 4, ## To show the numbers
	#'borderclrs'=> 'black',
	#do not draw border around graph
	#width of lines
	#'accent_treshold'=> 4,
	'legend_spacing' => 5, #spacing between legend elements
	#'legend_placement' =>'RC', #put legend to the centre right of chart
	#'dclrs' => \@catColor, #reference to array of category colours for each line
	'dclrs' => [ qw( green lred ) ], ## Fix the color ...
	'cumulate' => 'true', 	
	'transparent' => 0
	
)|| die "\nFailed to create cumulative bar chart: $my_graph->error()"; 
	

#set legend
$my_graph->set_legend(@allSps);
$my_graph->set_title_font(['verdana', 'arial', gdMediumBoldFont],16); 
$my_graph->set_x_label_font(['verdana', 'arial', gdMediumBoldFont],12); 
$my_graph->set_y_label_font(['verdana', 'arial', gdMediumBoldFont],12); 
$my_graph->set_x_axis_font(['verdana', 'arial', gdMediumBoldFont],12); 
$my_graph->set_y_axis_font(['verdana', 'arial', gdMediumBoldFont],12); 
$my_graph->set_legend_font(['verdana', 'arial', gdMediumBoldFont],12); 

#plot graph with table data
my $plot = $my_graph->plot(\@finalData);

open(IMG, ">$OutPutFile") || die ("\nFailed to save graph to file: $OutPutFile. $!");
print IMG $plot->gif();
close (IMG);
print "Created a cumulative bar chart for $resName Kb EBRs dataset in: $OutPutFile\n";

q^ if 0;
}

## Subroutines here --------------------------------------

sub readData {
my ($file, $file2)= @_;
my @allResolutions; my @allValues; my @allSpecies; my @decision;
open(FILE, "$file") || (warn "Can't open file $_\n");
      while (<FILE>) { 
	chomp;
	my @tmp=split(/\t/, $_); 
	push (@allValues, $_);
	push (@allSpecies, $tmp[2]);
	push (@decision, $tmp[5]);	
	}	
close(FILE);

my @uniqSps=uniq(@allSpecies);
my @decisionUniq=uniq(@decision);
my @uniqSpecies = sort { $a cmp $b } @uniqSps;
my @finalData; my @catColor; my @scores; my @all; my @uniqSpecies2; my $flag=0;

#foreach my $aaa(@uniqSpecies) { $aaa=~ s/_/ /ig; push (@uniqSpecies2, $aaa);}

#push (@finalData, \@uniqSpecies);
	foreach my $decision (@decisionUniq){
	foreach my $sps(@uniqSpecies) { 
		my $counter=0;
		foreach my $val (@allValues) { 
			my @line=split(/\t/, $val);	
			if (($sps eq $line[2]) and ($line[5] eq $decision)) { $counter++; }
		}
	if($counter == 0) { $counter=undef;}
	push (@scores, $counter); 
	$counter=0; my $newsps=$sps;
	$newsps=~ s/_/ /ig;
	push (@uniqSpecies2, $newsps);
	}
	if ($decision eq 'Break') { $decision = 'EBRs';} elsif ($decision eq 'PseudoBreak') { $decision = 'Gaps';} else { print '????\n';}
	push (@scores, $decision); ## To add the name of species
	my ($rand, $rand2);
	LABEL: { ($rand, $rand2) = random_colors(); redo LABEL if ($rand2 eq "white"); }
	if ($flag==0) { push(@finalData, [@uniqSpecies2]); $flag=1;} 
	push (@catColor, $rand); #print @catColor;
	push(@finalData, [@scores]);	
	undef @scores;
}
print2d (\@finalData, $file2);
return (\@finalData,\@catColor, \@decisionUniq); ## provided the breaks and pseudobreaks
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

