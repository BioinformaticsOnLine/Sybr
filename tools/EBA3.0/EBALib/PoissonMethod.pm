# $Id$ Poisson Method
# Perl module for EBA EBALib::PoissonMethod;
# Author: Jitendra Narayan <jnarayan81@gmail.com>
# Copyright (c) 2015 by Jitendra. All rights reserved.
# You may distribute this module under the same terms as Perl itself

##-------------------------------------------------------------------------##
## POD documentation - main docs before the code
##-------------------------------------------------------------------------##

=head1 NAME

EBALib::PoissonMethod  - DESCRIPTION of Object

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

package EBALib::PoissonMethod;

use strict;
#use warnings;
use Math::Round;
#use Term::ANSIColor;

use Exporter;

our @EXPORT_OK = "generateFinalPoisson";

sub generateFinalPoisson {
my ($path, $number, $lineage)= @_;
my $dir="$path/EBA_OutFiles";
EBALib::Messages::poissonMSG();

		 opendir(DIR, $dir) or die $!;
   		 while (my $file = readdir(DIR)) {
        		# We only want files
        		next unless (-f "$dir/$file");
        		# Use a regular expression to find files ending in .txt
        		next unless ($file =~ m/\.eba0$/); 
			my @new_file=split(/_/, $file);
				if ($new_file[0] ne "") { 
					my @res = split(/\//, $path);   EBALib::Messages::urRes($res[-1]); 
					my $InFileName="$path/EBA_OutFiles/final.eba6";
					my $AllFileName="$path/EBA_OutFiles/all_brk.eba0";
					my $OutFileName="$path/EBA_OutFiles/final_classify.eba7";
					poissonScore($InFileName, $AllFileName, $OutFileName, $res[-1], $number, $lineage);
				}				
		  }
closedir(DIR);
}

sub poissonScore {

my ($InFile, $InFile2, $OutFile, $resolution, $spsNumber, $lineage)=@_;

use strict;
#use warnings;
use Math::Round;

open(BETAFILE, 'betaScore') || (warn EBALib::Messages::failOp("betafile"));
my %betahash;
while (<BETAFILE>) { chomp;  my($key,$val) = split /\t/, lc($_); $_=EBALib::CommonSubs::trim($_); $key =~ s/^\s*(.*)\s*$/$1/; $val =~ s/^\s*(.*)\s*$/$1/; $betahash{$key} = $val; }
close BETAFILE or die EBALib::Messages::failCl("betafile");

$|++;
$/ = "\n";

open INFILE,  $InFile or die EBALib::Messages::failOp($InFile);
open OUTFILE, ">" , $OutFile or die EBALib::Messages::failOp($OutFile);

my @tmp_final; my %headhash; my $aa=0;  my @all_rj_scores; my @all_species; my @diff_species; my @fscore; my @all_lscore; my $gap_num=0; my $break_num=0; ### number of species is variable so we need to change

while (<INFILE>) {
   chomp;   
   @tmp_final = split /\t/, $_;
   s{^\s+|\s+$}{}g foreach @tmp_final;
   my $line= $_;

   #gap counting ===   temporary
   foreach my $vall(@tmp_final) { my @aa= split /\=/, $vall; my @bb= split /\+/, $aa[1];  if ($bb[0] eq "Gap"){ $gap_num++;} elsif ($bb[0] eq "Breakpoints"){ $break_num++;}}

   if ($. == 1 ) { foreach (@tmp_final) { 
	if ($aa <= ($spsNumber-1)) {$headhash{$_} = $aa; $aa++;}} 
	#print %headhash;  
	@all_species= keys %headhash; print OUTFILE "$line\n"; next;	
   }  ## So the header name should must be there, otherwise it will not able to find the species index.
   
   		open(CLASSFILE, "classification.eba") || (warn EBALib::Messages::failOp("classification.eba"));		
		while (<CLASSFILE>) {
   			chomp; my @all_sps_brk_cor;  my $Pvalue;  my @val2;  my @all_beta_score;   my $beta_score;  my @new_tmp_final; my $flag="OFF";
			if (index($_,"#") == 0) { next; } # Lines starting with a hash mark are comments
			if ($_ =~ /^\s*$/) { next; }
			$_=EBALib::CommonSubs::trim($_);
			my @tmp= split /\=/, lc($_);  
			next if $tmp[0] eq "classification"; 
			
			if ($tmp[0] eq "lineage") { @all_lscore = LineageScore (\@tmp_final, \%headhash, \%betahash, $resolution, $InFile2, $spsNumber); next;}
			my @class_sps = split /\,/, $tmp[1];    #print "$tmp[0]\n";
   			#print "@class_sps\n";
				foreach my $sps_name(@class_sps) {
					my $sps_index=$headhash{$sps_name};
					# print "$sps_name\t$sps_index\t$tmp_final[$sps_index]\n";
					my $sps_brk_cor=EBALib::CommonSubs::trim($tmp_final[$sps_index]);

					push @all_sps_brk_cor, $sps_brk_cor;
				
					my @sps_brk_cor_all = split /\=/, $sps_brk_cor;
					my @cor_info = split /\+/, $sps_brk_cor_all[1]; 
					next if $cor_info[0] eq "Gap";

						if ($sps_brk_cor eq 0) { 

							$beta_score = $betahash {"$resolution:$sps_name"};
							if($beta_score != 0) {push @all_beta_score, $beta_score;} ## I modified it .. becuase if there is any zero beta .. then it convert all to zero after multiplication. 
							#print "$beta_score======\n";
						} 
					 
						else    { 
							$beta_score = 1- $betahash {"$resolution:$sps_name"}; 
							push @all_beta_score, $beta_score;
							#print "$beta_score -----\n";
				        	}	
				}
			my $narrow_brk = getNarrowest (@all_sps_brk_cor);
			my $final_beta_score= EBALib::CommonSubs::multiply_all(@all_beta_score);
			undef @all_beta_score; undef @all_sps_brk_cor;
		        ## print "$tmp[0]\t$final_beta_score =====================>>>\n";
			# print "$line\t$narrow_brk\n";
			my %in_sps = map {$_ => 1} @class_sps;
			@diff_species  = grep {not $in_sps{$_}} @all_species;    ## to extract only those species whose name is not in the classification list. !!!
			# print "@diff_species\n";

				while (my ($key, $value) = each(%headhash)){
					my $brk_values=$tmp_final[$value];      ######hmmmmm @tmp_final??????????????????
				
					my @spsbrk= split /\=/, $brk_values;
					my @corinfo = split /\+/, $spsbrk[1]; 
					next if $corinfo[0] eq "Gap";
					# print "$corinfo[0]\n";

						if (grep /$key/, @diff_species) {
					 		if ($brk_values ne 0) {
	  
								my @breaks_cor = split /\=/, $brk_values;
								my @sp_cor = split /\,/, $breaks_cor[0];  @sp_cor = grep {$_} @sp_cor; ## it remove the 0 !!!
								# print "$sp_cor[0]\t$sp_cor[1]-----------------------------------------\n";
								my $nu=$spsNumber+1;   ## chromosome access from final file
								$Pvalue=getPoissonRate($key,$tmp_final[$nu],\@sp_cor, $InFile2,$spsNumber);   #print "$Pvalue\n";  .. 
                             					#print "$Pvalue\t$key\t$tmp_final[$nu]\t@sp_cor\t+++++++++\n";   
								my @val = split /\t/, $Pvalue;
								my @val2 = split /\t/, $narrow_brk;
								my $narrow_size=$val2[1]-$val2[0];
								my $result=$val[0]*($narrow_size + $val[1]);   ## print "$result=$val[0]*($narrow_size + $val[1])\n";    ### some of the values were in negative ... I need to check it !!!
								push @all_rj_scores, $result;
								$flag="ON";
							}
						}
				}
			
			my $final_rj_score= EBALib::CommonSubs::multiply_all(@all_rj_scores);
			undef @all_rj_scores;
			if ($flag eq "OFF") { $final_rj_score=1;}   #### need to improve !!!  In other word ... if breakpoints only in one species.
#print "$final_rj_score\t$flag---------------------------------->>>\n";
	                my $final = $final_rj_score * $final_beta_score;
##print OUTFILE "$final\t$final_rj_score\t$final_beta_score\t!!!!!!!!!!!\n";			
			push @fscore, "$tmp[0]:$final";
		undef @diff_species; undef %in_sps;
		} ##classification loop ends here 

		close CLASSFILE or die EBALib::Messages::failCl("file");

# my $firststring = join("\t",@fscore);
# my $firststring2 = join("\t",@all_lscore);

# print "$line\t$firststring\t$firststring2\n";

  push (@fscore, @all_lscore);
  s{^\s+|\s+$}{}g foreach @fscore;   # I remove leading and trailing whitespaces from array elements !!

  my %brkhash; my @brkscore; 

  if (($break_num == 1) and ($lineage)) { undef @fscore; $gap_num=0; $break_num=0;  next;}  ## Condition to not print lineage ... id "no" provided !!!!! ... now 0 values

  foreach my $brk_score (@fscore) { @brkscore = split /\:/, $brk_score; $brkhash{$brkscore[0]} = $brkscore[1]; }
  print OUTFILE "$line\t";
  my $count; my $first_val; my $second_val; my $ratio;
  print OUTFILE "@fscore\t";
	foreach my $value (sort {$brkhash{$b} <=> $brkhash{$a} } keys %brkhash) {    
           	$count++;
       		if ($count == 1) { print OUTFILE "$value:$brkhash{$value}\t"; $first_val=$brkhash{$value};}
       		if ($count == 2) { print OUTFILE "$value:$brkhash{$value}\t"; $second_val=$brkhash{$value};}
		if ($second_val !=0 ) { $ratio=$first_val/$second_val;} else { $ratio="NA";}
	}
  my $gap_brk_ratio= $gap_num/$spsNumber;  ### only for testing
  print OUTFILE "$ratio\t$gap_brk_ratio\t$spsNumber\t$gap_num\t$break_num\n";

 undef @fscore; $gap_num=0; $break_num=0; undef %brkhash; undef @brkscore; 
 #undef %betahash, undef %headhash;

} ## <INLFILE> loops ends here  

close INFILE or die EBALib::Messages::failCl("infile");
close OUTFILE or die EBALib::Messages::failCl("outfile");
 
} ## End main subroutine   !!!!!!!!   



## all other subroutines here

##===================================================================================================
## histogram calculation
sub histogram {
   my ($list_ref, $bin_width) = @_;
   my @list = @$list_ref;
   # This calculates the frequencies for all available bins in the data set
   my %histogram;
	
	$histogram{ceil(($_ + 1) / $bin_width) -1}++ for @list;
	
	return %histogram; 
undef %histogram;
}
##=====================================================================================================
# generating histogram 2d table

sub hist_table {
   my %histogram = %{shift()};
   my @hash_array;

   #use List::Util qw(max);
   # my $max = max values %histogram;

   my $max;
   my $min;
   # Calculate min and max
   while ( my ($key, $value) = each(%histogram) ) {
     $max = $key if !defined($min) || $key > $max;
     $min = $key if !defined($min) || $key < $min;
   }
	for (my $i = 0; $i <= $max; $i++) {
	my $frequency = $histogram{$i} || 0;
	$hash_array[$i]=$frequency;
	}
return @hash_array; 
undef @hash_array;
}
 
###------------------------------------------------------------------------------------------------------------------------
sub LineageScore {
my ( $tmp_final_ref, $headhash_ref, $betahash_ref, $resolution, $InFile2, $spsNumber) = @_;
    my %headhash = %$headhash_ref;
    my %betahash = %$betahash_ref;
    my @tmp_final = @$tmp_final_ref;
    my $key, my $value;		         
    my @all_species= keys %headhash; my $beta_sps;  my @all_scr; my $brk_size; my @all_final_res;

             while (($key, $value) = each(%headhash)){
		my $species_cor=$tmp_final[$value]; $species_cor =~ s/^\s*(.*)\s*$/$1/; 
		next if $species_cor == 0; ## To next if no breakpoints
		my @brk_cor = split /\=/, $species_cor;
                my @corinfo = split /\+/, $brk_cor[1];  
		next if $corinfo[0] eq "Gap";  ## Do not calculate the score for GAP !!!!!
		# my @sp_cor = split /\,/, $brk_cor[0];  @sp_cor = grep {$_} @sp_cor; ## it remove the 0 !!! Now no 0 in our dataset
		my @scord = split /\-\-/, $brk_cor[0];   ## what is more than one values !!!

		$brk_size= $scord[1]-$scord[0];
                       if ($species_cor ne 0) {				
				$beta_sps= 1-$betahash {"$resolution:$key"};
				my @diff_species = grep { $_ ne $key } @all_species;   # get the remaining !!!
					foreach my $sps_name (@diff_species) {
						my $sps_indx = $headhash {$sps_name};
			                        my $sps_cor=$tmp_final[$sps_indx]; 
				my @brk_cordi = split /\=/, $sps_cor; 
				my @cor_info = split /\+/, $brk_cordi[1]; 
				next if $cor_info[0] eq "Gap";
						if ($sps_cor ne 0) {
							my @breaks_cor = split /\=/, $sps_cor;
							my @sp_cor = split /\,/, $breaks_cor[0];  @sp_cor = grep {$_} @sp_cor; ## it remove the 0 !!!
							my $nu=$spsNumber+1;   ## chromosome access from final file
     							my $Pvalue=getPoissonRate($sps_name, $tmp_final[$nu],\@sp_cor, $InFile2, $spsNumber);   #print "$Pvalue\n"; 
                               				#print "$Pvalue=getPoissonRate($key, $tmp_final[$nu],@sp_cor,$InFile2, $spsNumber);\n"; 
							my @val = split /\t/, $Pvalue; 
							my $result=$val[0]*($brk_size + $val[1]);
							#print "$sps_name\t$result\t$val[0]\t($brk_size + $val[1])\n"; 
							push @all_scr, $result; 
		                                        }
			                        }
my $final_scr= EBALib::CommonSubs::multiply_all(@all_scr);
#print "$final_scr\t$final_scr\t$beta_sps\t===========================\n";
if (scalar (@all_scr) != 0 ) { $final_scr=$final_scr * $beta_sps;} else { $final_scr=$beta_sps; }  ## it assign zero if breaks only in one species.
#print "$key\t$final_scr\t==================@all_scr\n";
undef @all_scr;
push @all_final_res, "$key:$final_scr";
   			}
		}
#print "@all_final_res\n";
return @all_final_res;
undef @all_final_res; 
}   ## LineageScore subroutine ends here 

##-------------------------------------------------------------------------------------------------------------
# get the narrowest breakpoint amongst all
sub getNarrowest  {
    my @all_cor = @_; my @st_cor; my @ed_cor;
    foreach my $brk_cor (@all_cor) {
    	my @cor_val = split /\=/, $brk_cor;
	my @cor = split /\,/, $cor_val[0];
	foreach my $c(@cor) {
		next if $c == 0;
		my @cordi = split /\-\-/, $c;
		push @st_cor, $cordi[0]; 
		push @ed_cor, $cordi[1];
		}
	}

my @st_cor_sorted = sort { $a <=> $b } @st_cor;
my @ed_cor_sorted = sort { $a <=> $b } @ed_cor;
return "$st_cor_sorted[-1]\t$ed_cor_sorted[0]";
undef @st_cor, undef @ed_cor;
} 

# calculate the poisson rates and fectch them
sub getPoissonRate {
my ($sps_name, $chromo, $sps_cordi_ref, $InFile2, $spsNumber)= @_;
my @sps_cordi=@$sps_cordi_ref;

#print %hash;

use POSIX qw(ceil floor);
use List::Util qw(sum);

my @thefiles;  my @array_2d;   my @all_size;  my @sizes;  my @all_border; my @borders; my @all_chr;  my @all_chr2; my @chromosome; 
my $sizes_brk; my @sizes_brk;  my %hist_val;  my @hist_freq;  my @all_name; my $aa=0; my $chr; my $last_size; my @break_rates;

my @chr;
open(CHRFILE1, "chr_size.txt") || warn EBALib::Messages::failOp("chr_size");
while (<CHRFILE1>) { chomp; $_=EBALib::CommonSubs::trim($_); next if $_ =~ /^\s*#/; my @chrtmp = split /\t/, lc($_); push @chr,$chrtmp[0];}       ## We need to improve it !!!!!!!!!1
close CHRFILE1 or die EBALib::Messages::failCl("file");
for (@chr) { s/^\s+//; s/\s+$//;} # replace spaces
# my @chr=(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,"X");    ## Need to accept from the user command line or file. 

foreach $chr(@chr){
next if $chr ne $chromo;
        open INFILE2,  $InFile2 or die EBALib::Messages::failOp($InFile2);
	while (<INFILE2>) {                                             ### I think we can store the information in hash and used it  much faster than current !!!
	my $line=$_; chomp $line;
	my @tmp=split /\t/, $line;
	next if $sps_name ne $tmp[2];     ## Only read which match i.e file ...
	next if $tmp[5] ne "Break"; ## To ignore Pseudo breaks ..
	# my @name=split /\_/,$tmp[2];
		if (($chr eq $tmp[1]) and ($sps_name eq $tmp[2])) {
			my $size = $tmp[4]-$tmp[3];
			push (@all_size, $size);
			push (@all_border, "$tmp[3]\t$tmp[4]");
			push (@all_chr2, $tmp[1]); push (@all_chr, $tmp[1]);
			}
		} 
	close INFILE2 or die EBALib::Messages::failCl("file");
	
	my %histogram = histogram (\@all_size, 40000);	
	
	if ($last_size <= $#all_size) { $last_size=$#all_size};
        my @hist_table= hist_table(\%histogram);
	for (my $b=0; $b<=$last_size; $b++) { if($all_size[$b] eq "") { $all_size[$b]=0;} $sizes[$b][$aa]=$all_size[$b];}   

	for (my $dd=0; $dd<=$last_size; $dd++) { if($all_chr2[$dd] eq "") { $all_chr2[$dd]=0;} $chromosome[$dd][$aa]=$all_chr2[$dd];}

	for (my $cc=0; $cc<=$last_size; $cc++) { if($all_border[$cc] eq "") { $all_border[$cc]="0\t0";}  $borders[$cc][$aa]=$all_border[$cc]; }  
 
	for (my $bb=0; $bb<=$#hist_table; $bb++) { $hist_freq[$bb][$aa]="$chr:$hist_table[$bb]"; }
	undef @all_size; undef @all_border;  undef @all_chr2;  undef @all_chr; ### all undef here
$aa++;        	
}
my $chromosomeNumber=scalar(@chr);
my @break_rates2=fun_br_rates_chromo(\@hist_freq,\@sizes,$chromosomeNumber,40000, $spsNumber);   #number if chromosomes should be dynamic 

#print_2d (@break_rates2); ### Is to check the values after printing 
#print "$sps_name, $chromo, @sps_cordi\n";
my @scord = split /\--/, $sps_cordi[0];   ## what if more than one values !!! It is not possible to have more than one value .... !!!
my $brk_size= $scord[1]-$scord[0];
my @brk_rates2= map { $_->[0] } @break_rates2; ## convert to oneD !!!
 
#foreach (@brk_rates2) { print "$_\n";}
        my $indx=ceil($brk_size/40000); 
	my $PosValues= $brk_rates2[$indx-1];
	return "$PosValues\t$brk_size";
	
undef @chr; 
}   ## getPoissonRate subroutine ends here 


## list of Subrutines 

##===================================== CHROMOSOMES ==============================================
# Calculates breakpoint rates by bin non-homogeneous Poisson process;     
# Chromosme based R calculation .. I think so ... need to check in future !!!! 

sub fun_br_rates_chromo {
   my ($hist_br_ref,$size_br_ref,$Nchrom,$binsize, $spsNumber)= @_;
   my @hist_br = @$hist_br_ref;
   my @size_br = @$size_br_ref;
   
   my $values;  my @subArray;   my $number=0;  my @break_rates;   my $chr_size; my $jj=1;

my %hash; my @allChromo;
open(CHRFILE2, "chr_size.txt") || warn EBALib::Messages::failOp("chr_size");
while (<CHRFILE2>) { chomp; $_=EBALib::CommonSubs::trim($_); next if $_ =~ /^\s*#/; my ($key, $val) = split /\t/, lc($_); $hash{$key} = $val; push @allChromo, $key;}         ### We can read and store it ... !!!!
close CHRFILE2 or die EBALib::Messages::failCl("file");

my $maxGoes=scalar(@allChromo)*$spsNumber; undef @allChromo;

for (my $i=0; $i<=$#hist_br; $i++)
	{
	for (my $j=0; $j<=$maxGoes; $j++)  ## we need to change the number, as there the number of smepcies can varies.
		{   
			next if !$hist_br[$i][$j]; ####???????????????
			my $val=$hist_br[$i][$j]; next if $val eq ""; 
			my @both_values = split /:/, $val;  
			my @column = map { $$_[$j]} @size_br;
			foreach $values(@column) {
				my $ii= $i+1;
				if (($values> ($ii-1) * $binsize) && ($values < $ii * $binsize))
					{
					push (@subArray, $values);   #print "$values";
					}
				}
			 my $br_size_mean= mean (@subArray);
                         foreach (@column){ if ($_>0) { $number++;} }

		my $br_potential_sites = $hash{$both_values[0]} - 2*$br_size_mean/2; #exclude chromosome ends

		$br_potential_sites = $br_potential_sites - sum(@column) - $number*2*$br_size_mean/2; # exclude sites in and around other breakpoints 

                my $br_rates = $both_values[1]/$br_potential_sites;

$break_rates[$i][$j]=$br_rates;
undef @subArray; $number=0;
                      }   
	}
return @break_rates; 
undef %hash;
}            

# Calculate the mean 
sub mean { return @_ ? sum(@_) / @_ : 0 }

1;
