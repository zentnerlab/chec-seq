#!/usr/bin/perl

# Column	Value
# 0				Chromosome
#	1
# 2 			Start
# 3				End
# 4
# 5 			Length

die "perl pairs2bedgraph.pl <PAIRS FILE> <LOWER SIZE CUTOFF> <UPPER SIZE CUTOFF> <OUT FILE>\n" if(!$ARGV[3]);

print STDERR "Pairs file : $ARGV[0]\nSize Class : $ARGV[1] to $ARGV[2]\nOutput file : $ARGV[3].bedgraph\n";


$cutoff=500000;

open(FILE,$ARGV[0]) || die "INPUT $!\n";
while($line=<FILE>){
	@temp = split /[\ \s\n\t]+/, $line;
	#$temp[0]=~s/chr//;
	if($temp[5]>=$ARGV[1] && $temp[5]<=$ARGV[2]){
		for($i=$temp[2];$i<=$temp[3];$i++){
			$val{$temp[0]}{$i}++;
			$ln++;
		}
	}
	if( $ln >= $cutoff ){
		$cutoff+=500000;
		$str=sprintf("Reading %s\t%e\tCutoff %d\n",$temp[0],$ln,$cutoff) ;
		print STDERR $str
	}
}
close(FILE);

print STDERR "Done with pairs file\n";

#count total base-pairs mapped

foreach $i ( keys(%val) ){
	@x = keys(%{$val{$i}});
	$nbp += (1 + $#x);
}

print "$nbp\n";

open (OUT, ">$ARGV[3].bedgraph" ) || die "$!\n";

foreach $i ( keys(%val) ){
	foreach $j ( sort {$a <=> $b } ( keys(%{$val{$i}}) ) ){
		print OUT $i,"\t",$j,"\t",$j+1,"\t",($val{$i}{$j}*$nbp/$ln),"\n";
	}
}
