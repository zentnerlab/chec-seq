#!/usr/bin/perl -w

# Siva Kasinathan (Henikoff Lab, FHCRC): skasin@uw.edu

# Given a chromosome-specific BED file, a threshold value, and a max. distsance between
# distinct peaks, find all peak regions. Values are printed to STDOUT. Occupancies
# reported are areas under peaks.

# This script is an implementation of the general approach described in Kasinathan, et al.
# Nature Methods 2014.

use strict;
use autodie;

die "Usage: ./thresh_bed.pl [threshold] <path to chrom-specific bed> [interpeak distance]\n" if !$ARGV[2];

my $inter_peak = $ARGV[2];

my @vals; #This two-dimensional array will contain [coordinate,norm. count] pairs
my @vals_above_t; #This two-dimensional array will contain pairs that are above threshold

my $threshold = $ARGV[0];
my $chrom;

# Parse the chromosome-specific BED file
open BED, '<', $ARGV[1]
	or die "Could not open .gff file.";
	
while(<BED>){
	chomp;
	my @line = split /\t/,$_;
	push @vals,[$line[2],$line[3]];
	$chrom = $line[0];
}

close BED;

#Store positions above threshold in @vals_above_t
for (my $i = 0; $i < scalar @vals; $i++){
	if ($vals[$i][1] >= $threshold){
		push @vals_above_t,[$vals[$i][0],$vals[$i][1]];
	}
}

#Collapse positions that are nearby each other into peaks
my $counter = 0;
my $occupancy = 0;
my $start;
my @peaks; #This two-dimensional array will hold positions collapsed into peaks
#Each row of @peaks will contain three values (peak start, occupancy, peak width)
for (my $i = 0; $i < (scalar @vals_above_t)-1; $i++){
	$start = $vals_above_t[$i][0] if ($counter == 0);
	# If the next value is less than the inter peak distance away, consider it as
	# a part of the current peak
	if ($vals_above_t[$i+1][0] -$vals_above_t[$i][0] < $inter_peak){
		$counter++;
		$occupancy += $vals_above_t[$i][1];
	}
	# If the next value is greater than or equal to the inter peak distance away, consider
	# it a separate peak
	elsif ($vals_above_t[$i+1][0] -$vals_above_t[$i][0] >= $inter_peak){
		$counter++;
		$occupancy += $vals_above_t[$i][1];
		push @peaks,[$start,$occupancy,$vals_above_t[$i][0]];
		$occupancy = 0;
		$counter = 0;
	}
}

# Report all peak locations and occupancies
for (my $i = 0; $i < scalar @peaks; $i++){
	print $chrom, "\t", $peaks[$i][0], "\t", $peaks[$i][2]+1, "\t", $peaks[$i][1], "\n";
}

print STDERR "Num peaks: ", scalar @peaks, "\n";

exit;
