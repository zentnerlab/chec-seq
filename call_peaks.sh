#!/bin/bash

# Siva Kasinathan (Henikoff Lab, FHCRC): skasin@uw.edu

# This is a wrapper script that calls peaks using a thresholding approach
# Dependencies: BEDOPS (Neph, et al. Bioinformatics 2012; https://github.com/bedops/bedops)
#				treshold_bed.pl

# Usage: ./call_peaks.sh [threshold] [BED File] [Output prefix]

THRESH=$1	# Threshold
BEDFILE=$2	# Input BED file
OUTPREFIX=$3	# Output prefix

BEDOPSDIR=/home/gzentner/bedops
SCRIPTS_DIR=/home/gzentner/scripts

# Get list of chromosomes
CHROMS=(`${BEDOPSDIR}/bedextract --list-chr $BEDFILE`)

# Iterate over chromosomes
for c in ${CHROMS[@]}
do
	echo $c
	grep -w $c $BEDFILE > ${c}.temp	# Create temp file containing reads mapping to just a given chromosome
	for INTER in `seq 10 10 100`	# Iterate over a range of inter peak distances
	do
		# Call peaks using given parameters
		${SCRIPTS_DIR}/perl_thresh_bed.pl $THRESH ${c}.temp ${INTER} >> ${OUTPREFIX}_t${THRESH}_ip${INTER}.bedgraph
	done
done


rm -r *.temp	# Remove all temp files
