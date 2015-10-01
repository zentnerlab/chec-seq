#!/usr/bin/perl -w

# Siva Kasinathan (Henikoff Lab, FHCRC): skasin@uw.edu

# Score given sequences using a position-specific scoring matrix
# and filter based on p-value 

# Implementation of the algorithm underlying FIMO (Grant, Bailey, & Noble, Bioinformatics
# 2011) modified from code written for Univ. Wash. GENOME 541 course (taught by Bill Noble,
# Spring 2011).

# ScerTF motif format (tab delim; 4 rows x n columns where n is motif length)
# ScerTF motifs are available from: http://stormo.wustl.edu/ScerTF/download/
# ScerTF: Spivak & Stormo, NAR 2012.

use strict;
use autodie;
use POSIX;

die "Usage: pssm_scorer.pl <FASTA FILE> <MOTIF FILE> [p-value threshold]" if !$ARGV[2];

my $pval = $ARGV[2];
my @pssm;

#Read FASTA file into hash
my $seqs = &read_fasta($ARGV[0]);
my $num_seqs = scalar keys(%$seqs);
print STDERR "Successfully read FASTA file\n";

# Read position-specific scoring matrix and integerize
my($pssm,$motif_len) = &read_pssm($ARGV[1]);
print STDERR "Successfully read motif file\n";

#Implement dynamic programming to fill the score matrix and compute p-value array
my $p_vals = &calc_pvals($pssm,$motif_len);

my $strand_id;
my $max_score;
my $s_score;
my $rc_score;
my (@best_match,@substr_id,@strand,@best_score); #stores information for each sequence in FASTA File

for (my $i = 0; $i < $num_seqs; $i++){
	$best_match[$i] = $substr_id[$i] = $strand[$i] = $best_score[$i]= 0;
}

#Score each sequence and its reverse complement
my @ids = keys %$seqs;

for (my $i = 0; $i < $num_seqs; $i++){
	print STDERR "\t Processing sequence $i of $num_seqs\n" if ($i%1000 == 0);
	next if length($$seqs{$ids[$i]}) < $motif_len;
	my @substrings = kmers($$seqs{$ids[$i]},$motif_len);
	$max_score = -100000;
	for (my $k = 0; $k < $#substrings; $k++){
		my $s = $substrings[$k];
		#Store reverse complement of peak sequence
		my $rc_s = reverse $s;
		$rc_s =~ tr/ACGTacgt/TGCAtgca/;	
		
		#Score k-mer and its reverse complement
		$s_score = score($s,$pssm);
		$rc_score = score($rc_s,$pssm);
		
		#Figure out strand of motif match
		$strand_id = '+';	#Arbitrary assume best match on + strand, next check and fix this assumption if wrong
		if ($rc_score > $s_score){
			$s_score = $rc_score;
			$s = $rc_s;
			$strand_id = '-';
		}
		
		#Check if current k-mer or its rev. complement have better score than max score observed so far
		if ($s_score > $max_score){
			$max_score = $s_score;
			$best_match[$i] = $s;
			$substr_id[$i] = $k;
			$strand[$i] = $strand_id;
		}
	}
	$best_score[$i] = $max_score;
}


#Print output
#fasta_id	chr	motif_start	motif_end	matching_seq	strand	score	p_value
my @split_id;
my $below_thresh = 0;

print "Less than or equal to p-value threshold of ",$pval, "\n";
for (my $i = 0; $i < $num_seqs; $i++){
	if ($$p_vals[$best_score[$i]] <= $pval){
		$below_thresh++;
		@split_id = split_fasta_id($ids[$i]);
		$split_id[1] = 0 if (! defined $split_id[1]);
		print $ids[$i], "\t", $split_id[0], "\t", $split_id[1]+$substr_id[$i], "\t", $split_id[1]+$substr_id[$i]+$motif_len, "\t", $best_match[$i], "\t", $strand[$i], "\t", $best_score[$i], "\t", $$p_vals[$best_score[$i]],"\n" if $best_score[$i] > 0;
	}
}

# print "Greater than or equal to p-value threshold of ",$pval, "\n";
# for (my $i = 0; $i < $num_seqs; $i++){
# 	if ($$p_vals[$best_score[$i]] > $pval){
# 		@split_id = split_fasta_id($ids[$i]);
# 		$split_id[1] = 0 if (! defined $split_id[1]);
# 		print $ids[$i], "\t", $split_id[0], "\t", $split_id[1]+$substr_id[$i], "\t", $split_id[1]+$substr_id[$i]+$motif_len, "\t", $best_match[$i], "\t", $strand[$i], "\t", $best_score[$i], "\t", $$p_vals[$best_score[$i]],"\n" if $best_score[$i] > 0;
# 	}
# }

print "\n";
print "Percentage of sequences with p-value below threshold: ", $below_thresh, "/", $num_seqs, "=", $below_thresh/$num_seqs, "\n";
print "\n";

exit;

#Read FASTA file
sub read_fasta{
	open FASTA, "<", @_
		or die "Could not open FASTA file.\n";
	my %sequences;
	my $header;
	my $check = 0;
	my @lines = <FASTA>;
	
	foreach my $line (@lines){
		chomp($line);
		if($line =~ /^>/){
			$header = $line;
			$header =~ s/^>//;
			$header =~ s/\s.*//;
			if ($check == 0){
				$check = 1;
			}
			next;
		}
		if ($check == 0){
			die "File not in FASTA format.\n";
		}
		$line =~ tr/a-z/A-Z/;
		$sequences{$header}.=$line;
	}
	
	close(FASTA);
	return \%sequences;
} # End of sub read_fasta

#Read and integerize PSSM
sub read_pssm{
	open MATRIX, '<', @_
	or die "Could not open position-specific scoring matrix file.";

	my @pssm;

	while(<MATRIX>){
		chomp;
		my @line = split /\s+/,$_;
		shift @line;
		shift @line;
		push @pssm, [@line];
	}

	close MATRIX;
	
	#Find max and min entries of PSSM
	my $min = $pssm[0][0];
	my $max = $pssm[0][0];

	for(my $i=0; $i <= $#pssm; $i++){
		for(my $j=0; $j <= $#{$pssm[$i]}; $j++){
			$min = $pssm[$i][$j] if $pssm[$i][$j] < $min;
			$max = $pssm[$i][$j] if $pssm[$i][$j] > $max;
		}
	}

	#Calculate scaling factor
	$max=100/($max-$min);

	#Normalize by subtracting minimum value and multiplying by scaling factor
	for(my $i=0; $i <= $#pssm; $i++){
		for(my $j=0; $j <= $#{$pssm[$i]}; $j++){
			$pssm[$i][$j] -= $min;
			$pssm[$i][$j] *= $max;
			$pssm[$i][$j] = floor($pssm[$i][$j]+0.5);
		}
	}
	
	return (\@pssm,$#{$pssm[0]}+1);
} #End of sub read_pssm

#Calculate motif p-values
sub calc_pvals{
	my $pssm = $_[0];
	my $motif_len = $_[1];
	my @count_matrix;
	 
	for (my $yi=0; $yi < $motif_len; $yi++){
		for (my $xi=0; $xi <= 100*$motif_len; $xi++){
			$count_matrix[$yi][$xi] = 0;
		}
	}

	for (my $i=0; $i < 4; $i++){
		$count_matrix[0][${$pssm}[$i][0]]++;
	}

	for (my $i=1;$i < $motif_len; $i++){
		for (my $j=0; $j<=100*$motif_len; $j++){
			if ($count_matrix[$i-1][$j] != 0){
				for (my $k = 0; $k < 4; $k++){
					$count_matrix[$i][$j+$$pssm[$k][$i]] += $count_matrix[$i-1][$j];
				}
			}
		}
	}

	my $sum = 0;
	my @p_vals;

	for (my $i=0; $i<=100*$motif_len; $i++){
		$sum += $count_matrix[$motif_len-1][$i];
	}

	for (my $i=0; $i<100*$motif_len; $i++){
		push @p_vals,$count_matrix[$#count_matrix-1][$i]/$sum;
	}

	for (my $i=0; $i < $#p_vals; $i++){
		$sum = 0;
		for (my $j=$i; $j < $#p_vals; $j++){
			$sum += $p_vals[$j];
		}
		$p_vals[$i] = $sum;
	}
	
	return \@p_vals;
} #ENd of sub calc_pvals

#Calculate motif p-values
sub kmers{
	my($seq,$motif_length)=@_;	
	my @kseqs;
	
	for (my $i = 0; $i < length($seq) - $motif_length + 1; $i++){
		push @kseqs, substr $seq, $i, $motif_length;
	}
	
	return @kseqs;
}

#Score a k-mer given a pssm
sub score{
	my @bases=qw(A C G T);

	my($s,$pssm) = @_;
	$s =~ tr/a-z/A-Z/;
	my @seq = split //, $s;
	my $sum = 0;
	for (my $yi = 0; $yi < $#seq; $yi++){
		for (my $xi = 0; $xi < 4; $xi++){
			if ($seq[$yi] eq $bases[$xi]){
				$sum+=$$pssm[$xi][$yi];
				next;
			}
		}
	} 
	
	return $sum;
}

#Split FASTA ID formatted as chrXXX:start-end
sub split_fasta_id{
	my $id = $_[0];
	
	@split_id = split /[:,-]/, $id;
	
	return @split_id;
}
