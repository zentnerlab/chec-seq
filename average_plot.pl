#!/usr/bin/perl -w

# Siva Kasinathan (Henikoff Lab, FHCRC): skasin@uw.edu

# Given a regions BED file and BED file with positions falling in those regions with
# a particular score, return R/Excel-compatible file for calculation of statistics, e.g.,
# determination of aggregate profile

# The averages are printed to STDOUT while the raw data in R/Excel-compatible tabular form
# is printed to a designated output file.

use strict;
use autodie;

# Parse parameters

my ($regions_file, $bed_file, $output_file);
my $strand = 0;		#If strand orientation is present in regions BED file (fourth column),
					#then orient by strand

for (my $i = 0; $i < scalar @ARGV; $i++){
	if ($ARGV[$i] eq "-b"){
		$bed_file = $ARGV[$i+1];
		die "Bed file $bed_file does not exist!\n" unless (-e $bed_file);
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
	elsif ($ARGV[$i] eq "-h"){
		print_help();
	}
}

print_help() if (!defined($bed_file) || !defined($output_file) || !defined($regions_file));

# Read the BED file

my %bed;	# Hash indexed on chromosome and start coordinate that contains per-base counts
print STDERR "\tReading BED file\n";

my $n_bases = 0;	# Counter for number of bases

open BED, '<', $bed_file
	or die "Could not open BED file $bed_file\n";
while(<BED>){
	chomp;
	my @temp = split/\s+/;
	$bed{$temp[0]}{$temp[1]} = $temp[3];
	$n_bases++;
}
close BED;

print STDERR "Read $n_bases from $bed_file\n";


# Read in regions

my %region_bases;	# Hash indexed on name (chr_start_end_strand with strand being optional)
					# that contains the corresponding BED value read above if this value
					# exists
my $name;			# temp variable
my $n_regions = 0;	# counter for number of regions processed

print STDERR "\tProcessing regions .bed file\n";

open REGIONS, '<', $regions_file
	or die "Could not open regions BED file $regions_file\n";
while(<REGIONS>){
	chomp;
	my @temp = split/\s+/;
	$name = "$temp[0]"."_"."$temp[1]"."_"."$temp[2]";
	$name .= "_"."$temp[3]" if ($temp[3] eq "-" || $temp[3] eq "+");
	for (my $i = $temp[1]; $i < $temp[2]; $i++){
		$region_bases{$name}{$i} = $bed{$temp[0]}{$i} if (exists $bed{$temp[0]}{$i});
	}
	$n_regions++;
}
close REGIONS;

print STDERR "Read $n_regions from $regions_file\n";

# Output results

my @arr;	# Temporary array to hold all values in a window for a given locus/region

#Set field separator
local $" = "\t";

my (@sum,@n);	# Arrays to hold sum and number of values

print STDERR "\tOutputting results\n";

open OUTFILE, '>', $output_file
	or die "Could not open output file $output_file\n";
foreach my $region (sort keys %region_bases){
	my @temp = split(/_/,$region);	# Get name of region
	@arr=();	# reset array temp variable
	
	# Iterate through all bases in region
	for (my $i = $temp[1]; $i < $temp[2]; $i++){
		# If a value for this base in the region exists, store it in the array
		# if not, store "NA" in the array
		if (exists $region_bases{$region}{$i}){
			push @arr, $region_bases{$region}{$i};
		}
		else{
			push @arr, "NA";
		}
	}
	
	# Reverse values in the array if this locus/region is on the minus strand
	@arr = reverse @arr if ($strand && $temp[3] eq "-");
	
	# Calculate the sum and number of values for given positions in the window
	for (my $i = 0; $i < scalar @arr; $i++){
		if ($arr[$i] ne "NA"){
			$sum[$i] += $arr[$i];
			$n[$i]++;
		}
	}
	
	# Output values for each region (tabular format)
	print OUTFILE "@temp\t@arr\n";
	
}
close OUTFILE;

# Print averages

print STDERR "\tPrinting averages\n";

for (my $i = 0; $i < scalar @sum; $i++){
	if ($n[$i] > 0) {
		print $sum[$i]/$n[$i],"\n" 
	}
	else{
		print "NA\n";
	}
}

exit;

sub print_help{
	die "\nUsage: ./average_plot.pl [options] -b <bed_file> -r <regions_bed_file> -o <output_file>\n
	Required parameters:
	-b\tpath to .bed file
	-r\tpath to regions .bed file
	-o\tpath to output .txt file\n
	Optional parameters:
	-strand\torient by strand\n
	Help:
	-h Print usage\n\n";
}
