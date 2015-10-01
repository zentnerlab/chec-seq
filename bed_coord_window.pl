#!/usr/bin/perl -w

# Siva Kasinathan (Henikoff Lab, FHCRC): skasin@uw.edu

# Define windows around midpoints of BED coordinates. Note that if windows provided are
# not even integers, they are rounded up to the nearest even integer
# Output is to STDOUT

# This script performs an operation similar to SlopBED (in the BEDTools Suite)

use strict;
use autodie;

#Parse parameters
my ($bed_file, $output_file);
my $window;
my ($up,$down);

for (my $i = 0; $i < scalar @ARGV; $i++){
	if ($ARGV[$i] eq "-b"){
		$bed_file = $ARGV[$i+1];
		die "Bed file $bed_file does not exist!\n" unless (-e $bed_file);
	}
	elsif ($ARGV[$i] eq "-o"){
		$output_file = $ARGV[$i+1];
	}
	elsif ($ARGV[$i] eq "-w"){
		$window = $ARGV[$i+1];
		if (index($window,":") == -1){
			if (is_integer($window) && $window > 0){
				$window++ if $window % 2 != 0;
				$up = $down = $window/2;
			}
			else{
				die "\nSymmetric window width needs to be a positive integer\n";
			}
		}
		elsif (index($window,":") != -1){
			my @temp = split(/:/,$window);
			if (is_integer($temp[0]) && is_integer($temp[1]) && $temp[0] >= 0 && $temp[1] >= 0){
				$up = $temp[0];
				$down = $temp[1];
			}
			else{
				die "\nAsymmetric window width must be input as x:y where, for positive integers x and y,\nx is upstream 'half' window and y is downstream 'half' window\n\n";
			}
		}
	}
	elsif ($ARGV[$i] eq "-h"){
		print_help();
	}
}

print_help() if (!defined($bed_file) || !defined($output_file) || !defined($up) || !defined($down));

open BED, '<', $bed_file
	or die "Could not open .bed file $bed_file\n";
open OUTPUT, '>', $output_file
	or die "Could not open output file $output_file\n";

#Set field separator
local $" = "\t";

while(<BED>){
	chomp;
	my @temp = split/\s+/;
	$temp[1] += $temp[2];	# Add the start and end coordinates
	$temp[1]++ if $temp[1] % 2 != 0;	# Increment the sum of the start and end coordinates if not evenly divisible by 2
	$temp[1] /= 2;	# Calculate the midpoint of the feature
	$temp[2] = $temp[1]+1; # End of the midpoint is one base beyond the midpoint (BED convention)
	$temp[1] -= $up;	# Window (upstream)
	$temp[2] += $down;	# Window (downstream)
	print OUTPUT "@temp\n" if $temp[1] >= 0;
}


close BED;
close OUTPUT;

exit;

# Subroutine for checking if a string is an integer
sub is_integer {
   defined $_[0] && $_[0] =~ /^[+-]?\d+$/;
}

sub print_help{
	die "\nUsage: ./bed_coord_window.pl [options] -b <bed_file> -o <output_file> -w [window width]\n
	Outut is to STDOUT\n
	Required parameters:
	-b\tpath to .bed file
	-o\tpath to output BED file
	-w\twindow width (will be rounded if not an even integer)
	    [integer] will result in symmetric padding
	    [x:y] will result in asymmetric padding (x bp upstream, y bp downstream)\n
	Help:
	-h Print usage\n\n";
}
