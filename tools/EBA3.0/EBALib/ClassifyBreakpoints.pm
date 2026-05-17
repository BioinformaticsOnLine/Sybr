# $Id$ Classify Breakpoint
# Perl module for EBA EBALib::ClassifyBreakpoints;
# Author: Jitendra Narayan <jnarayan81@gmail.com>
# Copyright (c) 2015 by Jitendra. All rights reserved.
# You may distribute this module under the same terms as Perl itself

##-------------------------------------------------------------------------##
## POD documentation - main docs before the code
##-------------------------------------------------------------------------##

=head1 NAME

EBALib::ClassifyBreakpoints  - DESCRIPTION of Object

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

package EBALib::ClassifyBreakpoints;
use strict;
#use warnings;
#use Term::ANSIColor;

use Exporter;

our @EXPORT_OK = "classifyEBA";

sub classifyEBA {

my $reference=shift;
		
########################################
#
#  Perl Program to classify the species and convert them into evolutionary breakpoint format- Classify 
#  The location of taxdump folder ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
#
use strict;
#use warnings;


$|++; ## flush Perl print buffer

my (%names, %nodes, @array2d, %finalHash,@species, @allclass);

open SPSFILE, "sps.txt" or die $!;
open OUTFILE, ">", "classification.eba" or die $!;
while (<SPSFILE>) { chomp $_; @species=split /,/, lc($_); } 
close SPSFILE; ## It read the species names from sps.txt file ... need to improve !!!
unshift (@species, "$reference");
foreach (@species) { $_ =~ s/\s+/_/g;}   ## add underscore to space.

#######

my $taxonomy_directory='taxdump';
my ($files,$dirs)=getDirectoryFiles($taxonomy_directory);
process_file($_) for @$files;    # This subroutine create a hash of names and nodes file.

my %reverseNames;
while (my($key, $value) = each %names) { $reverseNames{lc($value)}=$key; }

my @inName=split(/:/, $reverseNames {lc($species[0])}); 
my $inputSpsId = $inName[0]; 
EBALib::Messages::refId($inputSpsId,$species[0]);

my $num=0; my %species_id;
EBALib::Messages::classMSG ();

foreach my $sp (@species) {
	my @array; 
	my @k=findkey ($sp, \%names );  
        s{^\s+|\s+$}{}g foreach @k; # Removing leading and trailing white space from array strings.
	my $lastsp=$sp; 
	# if (keys %{{ map {$_, 1} @k }} != 1) { @k[0]=senseMultipleHits(\@k,$sp);} 
	LINE: if (scalar(@k) > 1) { $k[0]=senseMultipleHits(\@k,$sp);}       
	elsif (!@k) { ($sp,@k)=senseNoHits(\@k,$sp); 
		foreach (0..$#species) { 
			if ($species[$_] eq "$lastsp") {s/$species[$_]/$sp/; }} 
		goto LINE;}
	else { print "\t@k\t$sp\n";}  # It search the "species names" in name file and extract its key.

 
	my @val = split(/:/, $k[0]);
	my $val = $val[0];  push @array, $val;  # Leaf of the tree ids are stored here, later we will add nodes in this array.
	$species_id{$val}=$sp;   # stored the hash for future use ... only studied species id 

		while ($val > 1) {
		my $parent = $nodes{$val};  # print "$parent\t";
        	#my $name= $names {$parent}; print "$name\n";  ### need to correct it print the different name bcz of duplication
		$val=$parent;
		push @array, $parent;
                }

	my @new_array=reverse @array;

	for (my $aa=0; $aa<= $#new_array; $aa++) { $array2d[$aa][$num]=$new_array[$aa]; }    ### Store the data in 2d array

undef @array; $num++;
}

%finalHash=classify(\@array2d, \%reverseNames);  ## the main subroutine to classify the species.

my @finalHashValues = values %finalHash;
my %finalHash2 = merge(\@finalHashValues, \%finalHash); 

my $counter=0;
foreach my $spsKey (keys %finalHash2) {       ## print the final result form finalHash array

    if ($counter == 0 ) { print OUTFILE "lineage=\n";}
    my $breakInSpecies=$finalHash2{$spsKey}; 
    my @breakInSpsArray= split /\,/ ,$breakInSpecies;
    my %in_brk = map {$_ => 1} @breakInSpsArray;
    my @diff  = grep {not $in_brk{$_}} @species;
if($spsKey=~ m/^[0-9]/) { my $vvv= $finalHash2{$spsKey}; print OUTFILE "$reference=$vvv\n"; next;}   ## Temporarily replacing the reference number with name ... need to improve later on !!!!
    print OUTFILE "$spsKey=$finalHash2{$spsKey}\n";
    push (@allclass, $spsKey);
   # @diff=join(",",@diff);  print "[!breaks=@diff]\n\n";
$counter++;
}
close OUTFILE or die EBALib::Messages::failCl("outfile");

# all subroutines here 

#----------------------------------------------
sub senseMultipleHits {
	my ($k_ref, $sp)=@_;
	my  @k=@$k_ref;
 	EBALib::Messages::multiHit();
	my %multipleId; 
	for(my $num=0; $num<=$#k; $num++) { print "\t$num\t$k[$num]\t$sp\n"; $multipleId{$num}=$k[$num];} 
	my $opted=<STDIN>; chomp($opted); 
	$k[0]=$multipleId{$opted}; 
	EBALib::Messages::opted($opted,$k[0],$sp);
return $k[0];
 ## Check the array for similarities of all hits and if they are different from one another then flash a message. 
}

#----------------------------------------------
sub senseNoHits {
  my ($k_ref, $sp)=@_;
  my  @k=@$k_ref;

  # sensorium to make sense of any errors
 EBALib::Messages::noHits($sp);
  my $newName = <STDIN>; chomp($newName); $newName=EBALib::CommonSubs::trim($newName); $newName =~ s/\s+/_/g;  ## replace the space with underscore ...
  my @spsId=findkey ($newName, \%names );

return ($newName,@spsId);
}

#----------------------------------------------------------------------
sub getDirectoryFiles{          # It get the directory files and return it
     my $taxdir = shift;

     opendir(my $dh, $taxdir) || die EBALib::Messages::failOp($taxdir);
     my @entries = grep {!( /^\.$/ || /^\.\.$/)} readdir($dh);
     @entries =  map { "$taxdir/$_" } @entries; #change to absolute paths
     closedir $dh;

     my @files =  grep( -f $_ , @entries);
     my @dirs = grep(-d $_, @entries);
     return (\@files,\@dirs);     ## return as a reference 

    close $dh;
}

#----------------------------------------------------------------------

sub process_file { # This is your custom subroutine to perform on each file    
    my $f = shift;     
    my ($val, $nam) = check_file($f);   
	if ($val == 1 and ($nam eq "names"))
		{  # print "processing file $f\n";
		%names=file2hash ($f, $nam);
		#return @array;
		}
         elsif ($val == 1 and $nam eq "nodes")
		{ # print "processing file $f\n";
		%nodes=file2hash ($f, $nam);
		}
        
	}

#-----------------------------------------------------------------------

sub check_file {
    use File::Basename;
    my $filepath = shift; # print $file;
    my $file = basename($filepath);
    my @ff= split /\./, $file;
    if ($ff[0] eq "names" || "nodes" ) 
	{ return 1, $ff[0]; }
}

#-------------------------------------------------------------------------
sub file2hash {
	my ($infile, $n) = @_;
	my %hash;  my %h;
	open FILE, $infile or die $!;
	while (<FILE>)
		{
   		chomp; my @tmpArray= split /\t\|\t/ , $_;   
		s{^\s+|\s+$}{}g foreach @tmpArray; # Removing leading and trailing white space from array strings.
		## next if $tmpArray[3] !~ m/^scientific/;
   		if ($n eq "names") { $tmpArray[0]="$tmpArray[0]:$.";}    # I make it unique by adding the line number and split later...
		$tmpArray[1] =~ s/\s+/_/g;  ## replace the space with underscore ...
		$hash{$tmpArray[0]} = $tmpArray[1]; 
	#	if($tmpArray[3] eq "scientific name") { $tmpArray[1] =~ s/\s+/_/g; $h{$tmpArray[0]} = $tmpArray[1];}  ## Make it global
		# print "$n\t$key\t$val\n";
        } 
	close FILE;
return %hash; 
     
}


#-------------------------------------------------------------------------
sub findkey {
        my ($species, $hash) =@_;
	my  %hash=%$hash;  my @all_keys;
	foreach my $key (keys %hash) {
     	if ($hash{$key} =~ m/^$species$/i) { push @all_keys, $key};
	}
s{^\s+|\s+$}{}g foreach @all_keys; # Removing leading and trailing white space from array strings.
return @all_keys;
undef @all_keys;
} 

#-------------------------------------------------------------------------
sub findId {
        my ($id, $hash) =@_;
	my  %hash=%$hash;  my @all_values;
	foreach my $key (keys %hash) {
	my @newValue= split(/:/, $hash{$key});       ## What if we have more than two hits for a key !!!!!
     	if ($newValue[0]==$id) { push @all_values, $key};
	}
s{^\s+|\s+$}{}g foreach @all_values; # Removing leading and trailing white space from array strings.
my $all_values=join(',',uniq(@all_values));
return $all_values;
undef @all_values;
} 

#---------------------------------------------------------------------
sub classify {
      my ($array2d_ref, $reverseNames_ref)=@_; 
      my  @array2d=@$array2d_ref;
      my  %reverseNames=%$reverseNames_ref;
      my @column;   my @all_leaf; my @row; my @uarray;  my @dup; my @all_g;  my $flag=0;

      for(my $i = 0; $i <= $#array2d; $i++){
	   for(my $j = 0; $j <= $#{$array2d[0]}; $j++){
	    push @row, $array2d[$i][$j];

		#for (@{$array2d[$i]}) { print "$_\t";}  ## it store all the values of $i rows.
		
		@uarray=uniq(@{$array2d[$i]}); 
		my %seen; @dup = map { 1==$seen{$_}++ ? $_ : () } @{$array2d[$i]};   ## To search only duplicated values in an array;
		@dup = grep {$_} @dup;   ## delete the blank or 0 in an array;
		# print "@dup\t---@uarray\t====\n";

        my @column = map {$$_[$j]} @array2d;   
		@column = grep {$_} @column; # delete the blank array at the end.
		push (@all_leaf, $column[-1]);
			}
            @all_leaf=uniq(@all_leaf);		
			#if (keys %{{ map {$_, 1} @dup }} == 1)
			if(scalar(@dup) == scalar(@uarray))
				{
				my @breaks_in= class2eba (\@all_leaf, \@all_leaf, $inputSpsId); 
				# print "@dup: @breaks_in\n";
				my @breakInNames =  returnName (\%species_id, \@breaks_in); # print @breakInNames;
				if ($flag != 1) {  my $brkin=join(",",uniq (@breakInNames)); $finalHash{$inputSpsId}= $brkin;}     
				$flag=1;
				}
 			else
				{
				my %all_clustered = clustered (\@dup,\@array2d);
				foreach my $group (keys %all_clustered) {
    				#print "The members of $group are\n";
					#no strict 'refs'; print "$group\n";
					next if !@{$all_clustered{$group}};
					my @uall=uniq(@{$all_clustered{$group}});
   					my @breaks_in= class2eba (\@uall, \@all_leaf, $inputSpsId);
					my @breakIn =  returnName (\%species_id, \@breaks_in);
					my $breaks_in=join(",",uniq (@breakIn)); ## print all sub breakpoints species
					my $groupName = findId ($group, \%reverseNames);   ### group is not identified in species_id because we have ids of only studied species.  Can use to extract name !!!!
					$finalHash{$groupName}= $breaks_in;
   					#print "$group:@breaks_in\n";
					push @all_g, $groupName;
    					}
				}
	undef @all_leaf; #print "$_\n";
	}

my $all_class = join(",",@all_g);
# $finalHash{'classification'}="lineage,$all_class,$inputSpsId";   # Not applicable in stand alone.
return %finalHash; # undef @all_g;
}
 #closedir(DIR);
##---------------------------------------------------------------------
sub merge {       # It merge the key from hash
        my ($values_array,$hash) =@_;
	my @values_array=@$values_array;
	my  %hash=%$hash;  
       	my %result; my @all_keys;
      
        for my $values (@values_array) {
		for my $key (keys %hash) {

        	my $accu = $hash{$key};
           	if ($values eq $accu) { push ( @all_keys, $key); }
	
            }
	my $keystring=join(',',@all_keys);

        undef @all_keys;
	$result{$keystring} = $values;
    }
return %result;
}

##------------------------------------------------------------------
sub returnName
{
my ($speciesId_ref,$breakIn_ref)= @_;
my %speciesId = %$speciesId_ref;
my @breakInName;
my @breakIn = @$breakIn_ref;
	foreach my $ids (@breakIn) { 
		my $idsName =  $speciesId { $ids };
		push (@breakInName, $idsName);
	#	print "$idsName\t";
		}
return @breakInName;
undef @breakInName;
}

#--------------------------------------------------------------------
	
sub uniq { my %seen; return grep { !$seen{$_}++ } @_; }


#-----------------------------------------------------------------------
sub clustered {
   my ($dup_ref,$array2d_ref)= @_;
   my @dup = @$dup_ref;
   my @array2d = @$array2d_ref;

   my @group; my %all_group;

	foreach my $d(@dup) { 

		for(my $i = 0; $i <= $#array2d; $i++){
	   		for(my $j = 0; $j <= $#{$array2d[0]} ; $j++){
	      			
               		 	my @column = map {$$_[$j]} @array2d;   
				@column = grep {$_} @column; # delete the blank array at the end.
				my $val = EBALib::CommonSubs::isInList ($d, @column);
				if ($val == 1) {
				push (@group, $column[-1]);  
						}
					}
				}
	$all_group{$d}=[@group]; 
	undef @group;			 
	}
return %all_group;
}

#----------------------------------------------------------------------------

sub class2eba  {
	my ($leaf_ref,$all_leaf_ref,$reference)= @_;
   	my @leaf = @$leaf_ref;
	my @all_leaf = @$all_leaf_ref;
	chomp($reference);
	my $val = EBALib::CommonSubs::isInList($reference, @leaf);
	 #print "$val\n";
		if($val != 1) { return @leaf;}
		else { 
	
			my %in_leaf = map {$_ => 1} @leaf;
			my @diff  = grep {not $in_leaf{$_}} @all_leaf;	
			       if ($#diff != -1) { return @diff; }
				else {
					my @new_leaf = grep { $_ != $reference } @all_leaf;
					return @new_leaf;
				     }
		    }

	}

#-----------------------------------------------------------------------------
# This function prints the script usage mode
sub printUsage {

  print << "End_Print_Usage";

Usage:
  perl ClassifyEBA.pl
------------------------------------------------------------------------------
Setting

   Need to download the taxdump zip folder from ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
   and unzip at the root. 

-------------------------------------------------------------------------------
Mandatory parameters:

  <Species Names> Name of the species should be comma sepeared. The of the species can be
   "Common Name" or "Scientific Name".

  <Reference Name> The first name of the input list will be treated as reference. 

  Warning:  Scientific name should be sepearted by only one space.

-------------------------------------------------------------------------------
Optional parameter:

  <RefName> If the first name is give as a star (*)  then it will cluster all the species in group.

-------------------------------------------------------------------------------
Output File formats:     Currently not generating !!!!

  + Breakpoint cassification outfile:

  The breakpoint classification outfile has the information about the breakpoints that were
  verified amongst genomes classofied according to reference species. 
	
  It has the following columns:     
  <Nodes> : Comma seperated species names
   File end with ";" sign.
--------------------------------------------------------------------------------
End_Print_Usage

exit(1);
}

} ## Main subroutines loop end here


1;
