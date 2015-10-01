#!/usr/bin/perl -w 

# Siva Kasinathan (Henikoff Lab, FHCRC): skasin@uw.edu

use strict;
use autodie;

# Given a set of 4-column BED files (i.e., bedgraph format – chromosome, start, end, value)
# Output a BED file wherein the maximum value across the input files is chosen for each position

die "Usage: ./combine_chec_bed.pl <BED file 1> <BED file 2> ...\nCan input any number of files OR a pattern\n" if scalar @ARGV==0;

#Hashes keyed on chromosome and base position that will contain
#Max value across datasets (in %merged) or file from which the max
#value was derived (in %source)... not the most memory-efficient implementation
my %merged;
my %source;

foreach my $file (@ARGV){ #Loop over all files passed to script through command line
	print STDERR $file,"\n";
	open BED, '<', $file
		or die "Could not open $file\n";
	
	while (<BED>){	#Parse each BED line
		chomp;
		my @line=split/\s+/;
		# If an entry does not exist in the hashes, create this entry
		if (!exists $merged{$line[0]}{$line[1]}){
			$merged{$line[0]}{$line[1]} = $line[3];
			$source{$line[0]}{$line[1]} = $file;
		}
		# If an entry exists in the hashes, update only if the value is greater than or equal
		# to the value that is contained in the hash
		elsif ($line[3] >= $merged{$line[0]}{$line[1]}){
			$merged{$line[0]}{$line[1]} = $line[3];
			$source{$line[0]}{$line[1]} .= ","."$file";
		}
	}
	
	close BED;
}

# Output results from merging in BED format
for my $chr (sort keys %merged){
	for my $coord (sort {$a <=> $b} keys %{$merged{$chr}}){
		#The following line is commented out to suppress printing of the source file for each value in the 5th colum of the output file.  This line can be uncommented if this information is desired.
		#print $chr, "\t", $coord, "\t", $coord+1, "\t", $merged{$chr}{$coord}, "\t", $source{$chr}{$coord},"\n";
		print $chr, "\t", $coord, "\t", $coord+1, "\t", $merged{$chr}{$coord}, "\n";
	}
}

exit;
