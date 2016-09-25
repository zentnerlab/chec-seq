#! /usr/bin/perl -w

# Siva Kasinathan (Henikoff Lab, FHCRC): skasin@uw.edu

# Take a pairs file (format below) and produce a BED file containing per-base counts
# of fragment ends. Optionally filter based on minimum and max. fragment lengths

# This script is adapted from a script by S. Ramachandran (Henikoff Lab, FHCRC)

# Column	Value
# 0			Chromosome
# 1			Strand
# 2 		Start
# 3			End
# 4			Map quality
# 5 		Length

die "./pairs_single_end_sizes.pl <PAIRS FILE> <OUT FILE> <LOWER BOUND> <UPPER BOUND>\n" if(!$ARGV[3]);

print STDERR "Pairs file : $ARGV[0]\nOutput file : $ARGV[1].bed\n";

#$cutoff=500000;

my $ln = 0;

open(FILE,$ARGV[0]) || die "INPUT $!\n";
while($line=<FILE>){
	@temp = split /[\ \s\n\t]+/, $line;
	
	# If a fragment is within the upper and lower bounds
	if ($temp[5] >= $ARGV[2] && $temp[5] <= $ARGV[3]){ 
		$val{$temp[0]}{$temp[2]}++;	# increment fragment start coordinate counts
		$val{$temp[0]}{$temp[3]}++;	# increment fragment end coordinate counts
		$ln+=2;	# increment number of mapped ends by two
	}
	
#	Print progress
# 	if( $ln >= $cutoff ){
# 		$cutoff+=500000;
# 		$str=sprintf("Reading %s\t%e\tCutoff %d\n",$temp[0],$ln,$cutoff) ;
# 		print STDERR $str;
# 	}
}
close(FILE);

print STDERR "Done with pairs file\n";

#count total base-pairs mapped

foreach $i ( keys(%val) ){
	@x = keys(%{$val{$i}});
	$nbp += (1 + $#x);
}

print STDERR "$nbp\n";

open (OUT, ">$ARGV[1].bed" ) || die "$!\n";

# Output normalized counts per base
# Normalization: multiple counts by number of bases mapped and divide by number of
# ends (i.e., number of fragments x 2)
foreach $i ( sort keys(%val) ){
	foreach $j ( sort {$a <=> $b } ( keys(%{$val{$i}}) ) ){
		print OUT $i,"\t",$j,"\t",$j+1,"\t",($val{$i}{$j}*$nbp/$ln),"\n";
	}
}

close OUT;

exit;
