#!/bin/bash

# Siva Kasinathan (Henikoff Lab, FHCRC): skasin@uw.edu

# Calculate averages over a set of genomic intervals with specified window width
# Dependencies: BEDOPS (Neph, et al. Bioinformatics 2012; https://github.com/bedops/bedops)
#				bed_coord_window.pl
#				average_plot.pl
#
# Note: This script is slightly different from window_ends.sh, which was written
#       specifically to deal with ends form paired-end sequencing data. This script
#       considers all bases contained within a paired-end read.

BEDOPS_DIR="/home/gzentner/bedops"
SCRIPTS_DIR="/home/gzentner/scripts"

display_usage(){
	echo ""
	echo "Usage: ./window.sh <BED file> <regions of interest> [window size]"
	echo ""
	echo "    <BED file>               path to BED file"
	echo "    <regions of interest>    regions of interest in BED format"
	echo "    [window size]            size of desired window in bp. Single positive integer will result in symmetric padding"
	echo "    <output file>            path to output file"
	echo ""	
}

if [ $# -ne 4 ]
then
	display_usage
	exit 1
fi

BED_FILE=$1
REGIONS=$2
WINDOW=$3
OUTFILE=$4

# Generate regions file with specified window parameters
echo "Determining windows around regions of interest"
WIN_REGIONS=$(basename "$REGIONS")
WIN_REGIONS=windowed.${WIN_REGIONS%.*}.tmp
${SCRIPTS_DIR}/bed_coord_window.pl -b $REGIONS -o $WIN_REGIONS -w $WINDOW

# Extract positions in BED file that fall within windowed regions
echo "Extracting BED positions that fall within windowed regions"
${BEDOPS_DIR}/bedextract $BED_FILE $WIN_REGIONS > temp

# Generate matrix
echo "Generating matrix for aggregate plots"
${SCRIPTS_DIR}/average_plot.pl -b temp -r $WIN_REGIONS -o ${OUTFILE} -strand
echo "Done";

rm temp
