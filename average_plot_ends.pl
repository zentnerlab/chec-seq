#!/usr/bin/perl -w

# Siva Kasinathan (Henikoff Lab, FHCRC): skasin@uw.edu

# Given a regions BED file and BED file with positions falling in those regions with
# a particular score, return R/Excel-compatible file for calculation of statistics, e.g.,
# determination of aggregate profile. This script specifically considers fragments ends
# that are either to the left of or to the right of a defined feature (strand information
# for the feature must be provided). The strand information of a read is not informative
# in this context.

# Averages are printed to STDOUT. R/Excel-compatible files containing arrays are printed
# to files

# Note: This script uses a 'pairs BED' file as input. This is a BED-formatted PAIRS file
# (essentially a reshuffling of columns of a PAIRS file).
# Column	Description
# 0			Chrom
# 1			Read start
# 2			Read end
# 3			Read length
# 4			Read quality
# 5			Strand

use strict;
use autodie;

# Parse parameters

my ($regions_file, $pairsbed_file, $output_file);
my $strand = 0;		#If strand orientation is present in regions BED file (fourth column),
					#then orient by strand

# Default read length minimum and maximum read lengths
my $min = 0;
my $max = 10000;

for (my $i = 0; $i < scalar @ARGV; $i++){
	if ($ARGV[$i] eq "-p"){
		$pairsbed_file = $ARGV[$i+1];
		die "Pairs.bed file $pairsbed_file does not exist!\n" unless (-e $pairsbed_file);
	}
	elsif ($ARGV[$i] eq "-o"){
		$output_file = $ARGV[$i+1];
	}
	elsif ($ARGV[$i] eq "-r"){
		$regions_file = $ARGV[$i+1];
		die "Regions .bed file $regions_file does not exist!\n" unless (-e $regions_file);
	}
	elsif ($ARGV[$i] eq "-strand"){
		$strand = 1;
	}
	elsif ($ARGV[$i] eq "-min"){
		$min = $ARGV[$i+1];
	}
	elsif ($ARGV[$i] eq "-max"){
		$max = $ARGV[$i+1];
	}
	elsif ($ARGV[$i] eq "-h"){
		print_help();
	}
}

print_help() if (!defined($pairsbed_file) || !defined($output_file) || !defined($regions_file));

# Read the input pairs BED file

# Hashes indexed on chromosome and start coordinate storing per-base fragments end counts
my %l_ends;
my %r_ends;
print STDERR "\tReading pairs.bed file\n";

open PAIRSBED, '<', $pairsbed_file
	or die "Could not open pairs.bed file $pairsbed_file\n";
while(<PAIRSBED>){
	chomp;
	my @temp = split/\s+/;
	# Store the fragment if its length is within the supplied min/max constraints
	if ($temp[3] >= $min && $temp[3] <= $max){
		$l_ends{$temp[0]}{$temp[1]}++;
		# Decrement the right ends by one (the end coordinate in a UCSC BED file is one
		# base beyond the actual feature end)
		$r_ends{$temp[0]}{$temp[2]-1}++;
	}
}
close PAIRSBED;

# Parse regions of interest

my %region_bases;	# Hash indexed on name (chrom_start_end_strand) with values
					# from left or right ends stored in an array (left is stored in pos 0
					# right is stored in pos 1)
# Temp variables:
my $name;
my ($right,$left);

print STDERR "\tProcessing regions .bed file\n";

open REGIONS, '<', $regions_file
	or die "Could not open regions BED file $regions_file\n";
while(<REGIONS>){
	chomp;
	# Decompose the region name
	my @temp = split/\s+/;
	$name = "$temp[0]"."_"."$temp[1]"."_"."$temp[2]";
	for (my $i = 3; $i <= $#temp; $i++){
		$name .= "_"."$temp[$i]" if ($temp[$i] ne "");
	}
	
	# For each position within the region
	for (my $i = $temp[1]; $i < $temp[2]; $i++){
		
		# Initialize the left and right read counts to a default value ("NA")
		$left = $right = "NA";
		
		# Get the left and right counts if they exist
		$left = $l_ends{$temp[0]}{$i} if exists $l_ends{$temp[0]}{$i};
		$right = $r_ends{$temp[0]}{$i} if exists $r_ends{$temp[0]}{$i};
		if ($temp[3] eq "-"){
			#If feature is on minus strand, then store right end as left end counts
			#and left end as right end counts
			push @{$region_bases{$name}{$i}}, $right;
			push @{$region_bases{$name}{$i}}, $left;
		}
		else{
			push @{$region_bases{$name}{$i}}, $left;
			push @{$region_bases{$name}{$i}}, $right;
		}
	}
}
close REGIONS;

# Temp arrays to contain left and right fragment end counts
my (@l_arr, @r_arr);
my ($rs,$ls);
#Set field separator
local $" = "\t";

# arrays containing sums and number of entries for the right and left reads
my (@l_sum,@l_n, @r_sum, @r_n);

print STDERR "\tOutputting results\n";

# Output files for left and right end arrays
my $l_output_file = "$output_file".".left.txt";
my $r_output_file = "$output_file".".right.txt";

# Output results
open LOUTFILE, '>', $l_output_file
	or die "Could not open output file $output_file\n";
open ROUTFILE, '>', $r_output_file
	or die "Could not open output file $output_file\n";
foreach my $region (sort keys %region_bases){	# Iterate over regions
	my @temp = split(/_/,$region);	# Split each region name
	# initialize empty temp arrays
	@l_arr=();
	@r_arr=();
	$rs = $ls = 0;
	
	# Iterate through each base in the region
	for (my $i = $temp[1]; $i < $temp[2]; $i++){
	# Save the values corresponding to those bases in the array
		push @l_arr, ${$region_bases{$region}{$i}}[0];
		push @r_arr, ${$region_bases{$region}{$i}}[1];
	}
	
	# If the feature is on the minus strand, reverse the arrays
	if ($temp[3] eq "-"){
		@l_arr = reverse @l_arr;
		@r_arr = reverse @r_arr;
	}
	
	# Compute sum and number of non-NA entries for left and right reads
	for (my $i = 0; $i < scalar @l_arr; $i++){
		if ($l_arr[$i] ne "NA"){
			$l_sum[$i] += $l_arr[$i];
			$l_n[$i]++;
			$ls++;
		}
		if ($r_arr[$i] ne "NA"){
			$r_sum[$i] += $r_arr[$i];
			$r_n[$i]++;
			$rs++;
		}
	}
	
	# Print the array values for a given region to the output file if there is at least
	# one non-NA value
	print LOUTFILE "@temp\t@l_arr\n" if $ls > 0;
	print ROUTFILE "@temp\t@r_arr\n" if $rs > 0;
	
}
close LOUTFILE;
close ROUTFILE;


#Print left and right averages to STDOUT
#Output format: [left_avg] [right_avg]
for (my $i = 0; $i < scalar @l_sum; $i++){
	if ($l_n[$i] > 0){
		print $l_sum[$i]/$l_n[$i];
	}
	else{
		print "NA";
	}
	print "\t";
	if ($r_n[$i] > 0){
		print -1*$r_sum[$i]/$r_n[$i];
	}
	else{
		print "NA";
	}
	print "\n";
}

exit;

sub print_help{
	die "\nUsage: ./average_plot_ends.pl [options] -b <bed_file> -r <regions_bed_file> -o <output_file>\n
	Required parameters:
	-p\tpath to .bed file
	-r\tpath to regions .bed file
	-o\toutput file PREFIX\n
	Optional parameters:
	-strand\torient by strand\n
	Help:
	-h Print usage\n\n";
}
