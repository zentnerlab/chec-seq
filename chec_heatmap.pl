#!/usr/bin/perl -w

# S. Kasinathan, Henikoff Lab/FHCRC (skasin@uw.edu)

# Get values corresponding to ChEC time points falling within defined regions
# (peaks) and output in format compatible with heatmap generation

# Input:  Bedgraph containing peaks, parameters file pointing to time point bedgraph files
# Output: Array containing peak, coordinates, occupancy and desired statistic for the peak
#         region across ChEC experiments/time points

# Notes on processing of input data in this script:
# (1) Z-score transformation: In order to facilitate comparisons between peaks that have
# different occupancies over time, values within a heatmap row (with each entry in a row
# corresponding to a ChEC timepoint) are Z-score transformed. Z-score transformation can
# be suppressed with the '-no_z' flag.
# (2) Background normalization: This corrects for signal:noise ratio and allows a more
# robust comparison between ChEC timepoints. Specifically, without this normalization, the 
# early ChEC timepoints tend to have high normalized counts at peaks due to low background
# cleavage, while late time points have comparatively lower normalized counts at peaks
# due to higher background. In order to facilitate visual comparison of timepoints in
# a heatmap without the potentially confounding effect of the particulars of signal:noise,
# normalization relative to the immediate local background around a peak is performed.
# Background window width can be modulated using the '-bg_w' flag (default width is 25 bp
# upstream and downstream of the peak of interest). Background normalization can be
# turned off with the '-no_bg' flag. Note: This feature was not utilized in the manuscript.

# Type "./chec_heatmap.pl" to display usage.

use strict;
use autodie;

print_help() if !$ARGV[0];

#-----------------------------------------------------------------------------------------
# Parse parameters
#-----------------------------------------------------------------------------------------

my $peaks_file;
my $file_list;

my $calc_z = 1;     # Whether or not to calculate Z-scores (default: calculate z-scores)
my $stat ="avg";    # Statistic to compute (default: average)
my $min_t_vals=3;	# Minimum number of time points for which there needs to be data
                    # default: 3

my $calc_bg = 1;	# Binary variable indicating whether or not perform background
                    # normalization. Default: yes
my $bg_width = 25;  # Default background size. 25 bp upstream of peak and 25 bp downstream

for (my $i = 0; $i < scalar @ARGV; $i++){
	if ($ARGV[$i] eq "-p"){
		$peaks_file = $ARGV[$i+1];
		die "Pairs file $peaks_file does not exist!\n" unless (-e $peaks_file);
	}
	elsif ($ARGV[$i] eq "-d"){
		$file_list = $ARGV[$i+1];
		die "Data file list $file_list does not exist!\n" unless (-e $file_list);
	}
	elsif ($ARGV[$i] eq "-s"){
		$stat = $ARGV[$i+1];
	}
	elsif ($ARGV[$i] eq "-no_z"){
		$calc_z = 0;
	}
	elsif ($ARGV[$i] eq "-min_t"){
		$min_t_vals = $ARGV[$i+1];
		die "\nMin. number of data points needs to be an integer\n\n" if (!is_integer($min_t_vals) || $min_t_vals <= 0);
	}
	elsif ($ARGV[$i] eq "-bg_w"){
		$bg_width = $ARGV[$i+1];
		die "\nBackground width needs to be positive integer\n\n" if (!is_integer($bg_width) || $bg_width <= 0);
	}
	elsif ($ARGV[$i] eq "-no_bg"){
		$calc_bg = 0;
	}
	elsif ($ARGV[$i] eq "-h"){
		print_help();
	}
}

print STDERR "\nParameters:\n";
print STDERR "\tChEC timepoint datasets file:   $file_list\n";
print STDERR "\tPeaks file:                     $peaks_file\n";
print STDERR "\tMinimum number of datapoints:   $min_t_vals\n";
print STDERR "\tStatistic:                      $stat\n";
print STDERR "\tCalculate z-score:              ";
if ($calc_z){
	print STDERR "yes\n";
}
else{
	print STDERR "no\n";
}
print STDERR "\tBackground normalization:       ";
if ($calc_bg){
	print STDERR "yes\n";
	print STDERR  "\tBackground width:               $bg_width\n";
}
else{
	print STDERR "no\n";
}
print STDERR "\n";


#-----------------------------------------------------------------------------------------
# (1) Read peaks bedgraph
#-----------------------------------------------------------------------------------------

print STDERR "Reading peaks bedgraph $peaks_file\n";

my $npeaks = 0; # Counter for number of peaks
my $nbases = 0; # Counter for number of bases contained within peaks

# Hashes keyed by chromosome and start position relating to peaks
my %peak_bases;	# Bases contained within peaks
my %peaks;		# Chrom and start and end coordinates of peaks
my %occupancies; # Occupancies of peaks (or other values provided in input bedgraph
                 # corresponding to each peak.

# Hashes keyed by chromosome and start position relating to background
my %bg_bases;

open PEAKS, '<', $peaks_file
	or die "Could not open peaks bedgraph $peaks_file\n";

while (<PEAKS>){
	chomp;
	my @line = split/\s+/;
	$peaks{$line[0]}{$line[1]} = $line[2];
	$occupancies{$line[0]}{$line[1]} = $line[3];
	for (my $pos = $line[1]; $pos < $line[2]; $pos++){
		$peak_bases{$line[0]}{$pos} = 1;
		$nbases++;
	}
	if ($calc_bg){
		for (my $pos = $line[1] - $bg_width; $pos < $line[1]; $pos++){
			$bg_bases{$line[0]}{$pos} = 1;
		}
		for (my $pos = $line[2]; $pos < $line[2] + $bg_width; $pos++){
			$bg_bases{$line[0]}{$pos} = 1;
		}
	}
	$npeaks++
}

close PEAKS;

print STDERR "\tDone. Read $npeaks peaks ($nbases total bases).\n";
print STDERR "\n";

#-----------------------------------------------------------------------------------------
# (2) Read time point bedgraph files
#-----------------------------------------------------------------------------------------

print STDERR "Parsing timepoint dataset list file and reading datasets\n";

my $nfiles = 0; # Counter for number of timepoint bedgraph files (ChEC cleavage datasets)
my %tfiles;

open FILE_LIST, '<', $file_list
	or die "Could not open list of ChEC cleavage datasets $file_list\n";

while(<FILE_LIST>){
	chomp;
	my @line = split/\s+/;
	$tfiles{$line[0]} = $line[1];
	$nfiles++;
}

close FILE_LIST;

print STDERR "\tDone. Total of $nfiles ChEC datasets to process\n";

my %data;
my $ntdata = 0; # Counter for number of bases from a given timepoint that fall within peaks
my $nbg_data = 0; # Counter for number of bases from a given timepoint that fall in
                  # peak background regions
                  
for my $time (sort {$a <=> $b} keys %tfiles){
	open TIMEPOINT, '<', $tfiles{$time}
		or die "Could not open ChEC dataset $tfiles{$time}\n";
	
	print STDERR "\tProcessing time point $time (file: $tfiles{$time})\n";
	$ntdata = 0;	
	while(<TIMEPOINT>){
		chomp;
		my @line = split/\s+/;
		# Check if data point is in a peak; if it is, update the hash containing
		# time point data.
		if (exists $peak_bases{$line[0]}{$line[1]}){
			$data{$line[0]}{$line[1]}{$time} = $line[3];
			$ntdata++;
		}
		if ($calc_bg && exists $bg_bases{$line[0]}{$line[1]}){
			$data{$line[0]}{$line[1]}{$time} = $line[3];
			$nbg_data++;
		}
	}
	
	close TIMEPOINT;
	print STDERR "\t\tDone. Number of bases from timepoint falling in peaks: $ntdata\n";
	print STDERR "\t\t      Number of bases from timepoint falling in background: $nbg_data\n" if $calc_bg;
}

print STDERR "\tDone.\n\n";

#-----------------------------------------------------------------------------------------
# (3) Compute statistic of interest for peaks through time
#-----------------------------------------------------------------------------------------

print STDERR "Computing $stat across time points for peaks and background.\n";

my %time_res;
my %sig_noise;

my %temp; # Temporary hash to hold values over all experiments for bases in a given peak
my $val; # Temp variable to hold statistic computed for a given time point
my @hash_to_arr;
for my $chr (sort keys %peaks){	# Loop over all chromosomes
	for my $start (sort {$a <=> $b} keys %{$peaks{$chr}}){ # Loop over all start positions
		my $end = $peaks{$chr}{$start};	# Get peak end coordinate
		%temp=(); # Clear temp hash
		for (my $pos = $start; $pos < $end; $pos++){	# Loop over all positions in peak
			if (exists $data{$chr}{$pos}){	# If there is data for this position in peak
				# Loop over all time point positions for which there is data:
				for my $tpoint (keys %{$data{$chr}{$pos}}){
					$temp{$tpoint}{$pos} = $data{$chr}{$pos}{$tpoint};
				}
			}
		}	#End peak positions loop
	
		# Compute statistic of interest
		for my $tpoint (sort {$a <=> $b} keys %temp){
			# Convert data for bases to ordered list (array)
			@hash_to_arr = values %{$temp{$tpoint}};

			# Compute desired statistic
			if ($stat eq "avg"){
				$val = mean(\@hash_to_arr);
			}
			elsif ($stat eq "sum"){
				$val = sum(\@hash_to_arr);
			}
			elsif ($stat eq "min"){
				$val = min_val(\@hash_to_arr);
			}
			elsif ($stat eq "max"){
				$val = max_val(\@hash_to_arr);
			}
			elsif ($stat eq "stdev"){
				$val = stdev(\@hash_to_arr);
			}
		
			$time_res{$chr}{$start}{$tpoint} = $val;
			$sig_noise{$tpoint}{$chr}{$start} = sum(\@hash_to_arr);
		}
	
	} # End peak start coordinates loop
} # End chromosomes loop

my %bg_norm;

if ($calc_bg){
	for my $chr (sort keys %peaks){	# Loop over all chromosomes
		for my $start (sort {$a <=> $b} keys %{$peaks{$chr}}){ # Loop over all start positions
			my $end = $peaks{$chr}{$start};	# Get peak end coordinate
			%temp=(); # Clear temp hash
			# Get left half of background
			for (my $pos = $start-$bg_width; $pos < $start; $pos++){	# Loop over all positions in peak
				if (exists $data{$chr}{$pos}){	# If there is data for this position in peak
					# Loop over all time point positions for which there is data:
					for my $tpoint (keys %{$data{$chr}{$pos}}){
						$temp{$tpoint}{$pos} = $data{$chr}{$pos}{$tpoint};
					}
				}
			}	#End peak positions loop
		
			#Get right half of background
			for (my $pos = $end; $pos < $end+$bg_width; $pos++){	# Loop over all positions in peak
				if (exists $data{$chr}{$pos}){	# If there is data for this position in peak
					# Loop over all time point positions for which there is data:
					for my $tpoint (keys %{$data{$chr}{$pos}}){
						$temp{$tpoint}{$pos} = $data{$chr}{$pos}{$tpoint};
					}
				}
			}	#End peak positions loop
	
			# Compute statistic of interest
			for my $tpoint (sort {$a <=> $b} keys %temp){
				# Convert data for bases to ordered list (array)
				@hash_to_arr = values %{$temp{$tpoint}};

				$sig_noise{$tpoint}{$chr}{$start} /= sum(\@hash_to_arr) if exists $sig_noise{$tpoint}{$chr}{$start};
			}
	
		} # End peak start coordinates loop
	} # End chromosomes loop

	# Compute average signal:noise
	for my $tpoint (keys %tfiles){
		@hash_to_arr=();
		for my $chr (keys %{$sig_noise{$tpoint}}){
			for my $start (keys %{$sig_noise{$tpoint}{$chr}}){
				push @hash_to_arr, $sig_noise{$tpoint}{$chr}{$start};
			}
		}
		$bg_norm{$tpoint} = mean(\@hash_to_arr);
		print STDERR $tpoint, "\t", $bg_norm{$tpoint}, "\n";
	}

	@hash_to_arr = ();
	@hash_to_arr = values %bg_norm;
	my $max = max_val(\@hash_to_arr);
	for my $tpoint (keys %bg_norm){
		$bg_norm{$tpoint} /= $max;
	}

#	Perform background normalization:
# 	for my $chr (keys %time_res){
# 		for my $start (keys %{$time_res{$chr}}){
# 			for my $tpoint (keys %{$time_res{$chr}{$start}}){
# 					if ($time_res{$chr}{$start}{$tpoint} ne "NA"){
# 						$time_res{$chr}{$start}{$tpoint} /= $bg_norm{$tpoint};
# 				}
# 			}
# 		}
# 	}

}

print STDERR "\tDone.\n\n";

#-----------------------------------------------------------------------------------------
# (4) Z-score transform time-resolved peak data
#-----------------------------------------------------------------------------------------
# Loop over hash containing statistic of interest for time points
if ($calc_z){

	print STDERR "Performing Z-score transformation and filtering data to meet minimum criterion.\n";
	my $excl = 0;
	my $incl = 0;
	
	my $ctr = 0;
	for my $chr (keys %time_res){
		for my $start (keys %{$time_res{$chr}}){
			# Convert data for time points to ordered list (array)
			@hash_to_arr =  @{$time_res{$chr}{$start}}{sort {$a <=> $b} keys %{$time_res{$chr}{$start}}};
		
			$ctr = 0;
			if (scalar @hash_to_arr >= $min_t_vals){
				$incl++;
				z_score(\@hash_to_arr);
				
				for my $tpoint (sort {$a <=> $b} keys %{$time_res{$chr}{$start}}){
					$time_res{$chr}{$start}{$tpoint} = $hash_to_arr[$ctr];
					$ctr++;
				}
			}
			else{
				delete $time_res{$chr}{$start};
				$excl++;
			}
		}
	}
	print STDERR "\tDone. Number of bases that did not have at least $min_t_vals datapoints: $excl.\n";
	print STDERR "\t      Number of bases that have at least $min_t_vals datapoints: $incl.\n\n";
}

#-----------------------------------------------------------------------------------------
# (5) Output results
#-----------------------------------------------------------------------------------------

print STDERR "Outputting results.\n";

for my $chr (sort keys %time_res){
	for my $start (sort {$a <=> $b} keys %{$time_res{$chr}}){
		print $chr, "\t", $start, "\t", $peaks{$chr}{$start}, "\t", $occupancies{$chr}{$start};
		for my $tpoint (sort {$a <=> $b} keys %tfiles){
			if (exists $time_res{$chr}{$start}{$tpoint}){
				print "\t", $time_res{$chr}{$start}{$tpoint};
			}
			else{
				print "\tNA";
			}
		}
		print "\n";
	}
}

print STDERR "\tDone.\n\n";

print STDERR "Job done.\n\n";

exit;

#-----------------------------------------------------------------------------------------
# Subroutines
#-----------------------------------------------------------------------------------------

sub print_help{
	print "\nUsage: ./chec_heatmap.pl [options] -p <peaks file> -d <ChEC bedgraphs>\n\n";
	print "-p     peaks file (bedgraph)\n";
	print "-d     two-column fine containing time (in seconds) and path to\n";
	print "       ChEC dataset in bedgraph format\n";
	print "\nOptions:\n";
	print "-s     statistic to compute; options: sum (default), max, min, mean, median\n";
	print "-bg_w  background width (e.g., \"-bg_w 50\" is 50 bp on either side of peak);\n";
	print "       default: 25 bp\n"; 
	print "-no_bg turn off background normalization (default: on)\n";
	print "-no_z  turn off z-score calculation (default: compute z-score)";
	print "-min_t minimum number of input ChEC experiments for which data needs to be present\n";
	print "       in order to be be outputted (default: 3)\n\n";
	print "Help:\n";
	print "-h     print this help screen\n\n";
	
	exit;
}

# Check if a value is an integer
sub is_integer {
   defined $_[0] && $_[0] =~ /^[+-]?\d+$/;
}

# Return the sum of an array passed by reference
sub sum{
	my $arr = shift;
	my $sum = 0;
	
	if (scalar @{$arr} == 0){
		$sum = "NA";
	}
	else{
		for (my $i = 0; $i < scalar @{$arr}; $i++){
			$sum += $$arr[$i] if ($$arr[$i] ne "NA");
		}
	}
	
	return $sum;
}

# Return the minimum value in an array passed by reference
sub min_val{
	my $arr = shift;
	my $min;
	
	if (scalar @{$arr} == 0){
		$min = "NA";
	}
	else{
		$min = $$arr[0];
		for (my $i = 0; $i < scalar @{$arr}; $i++){
			$min = $$arr[$i] if ($$arr[$i] < $min);
		}
	}
	
	return $min;
}

# Return the maximum value in an array passed by reference
sub max_val{
	my $arr = shift;
	my $max;
	
	if (scalar @{$arr} == 0){
		$max = "NA";
	}
	else{
		$max = $$arr[0];
		for (my $i = 0; $i < scalar @{$arr}; $i++){
			$max = $$arr[$i] if ($$arr[$i] < $max);
		}
	}
	
	return $max;
}

# Compute the mean of an array passed by reference
sub mean{
	my $arr = shift;
	my $mean;
	
	if (scalar @{$arr} == 0){
		$mean = "NA";
	}
	else{
		my $n = 0;
		for (my $i = 0; $i < scalar @{$arr}; $i++){
			$n++ if $$arr[$i] ne "NA";
		}
	
		$mean = sum($arr)/$n;
	}
	
	return $mean;
}

# Compute the standard deviation of an array passed by reference
sub stdev{
	my $arr = shift;
	my $sigma;
	
	if (scalar @{$arr} == 0){
		$sigma = "NA";
	}
	else{
		my $mu = mean($arr);
		my $n = 0;
			for (my $i = 0; $i < scalar @{$arr}; $i++){
				$n++ if $$arr[$i] ne "NA";
			}
	
		if ($n > 1){
			for (my $i = 0; $i < scalar @{$arr}; $i++){
				$sigma += ($$arr[$i]-$mu)**2 if $$arr[$i] ne "NA";
			}
			$sigma = sqrt($sigma/($n-1));
		}
		else{
			$sigma = "NA";
		}
	}
	
	return $sigma;
}

# Compute Z-score of an array passed by reference
sub z_score{
	my $arr = shift;
	my $sigma = stdev($arr);
	my $mu = mean($arr);
	
	for (my $i = 0; $i < scalar @{$arr}; $i++){
		if ($$arr[$i] ne "NA" && $sigma ne "NA"){
			$$arr[$i] = ($$arr[$i] - $mu)/$sigma;
		}
		else{
			$$arr[$i] = "NA";
		}
	}
	return;
}

# EOF #
