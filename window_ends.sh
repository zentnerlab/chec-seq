#!/bin/bash

# Siva Kasinathan (Henikoff Lab, FHCRC): skasin@uw.edu

# Calculate averages over a set of genomic intervals with specified window width
# Dependencies: BEDOPS (Neph, et al. Bioinformatics 2012; https://github.com/bedops/bedops)
#				bed_coord_window.pl
#				average_plot_ends.pl

# Directory containing BEDOPS executables				
BEDOPS_DIR="/home/gzentner/bedops"

# Directory containing scripts (bed_coord_window.pl and average_plot_ends.pl)
SCRIPTS_DIR="/home/gzentner/chec_paper_scripts"

display_usage(){
	echo ""
	echo "Usage: ./window_ends.sh -min [read len] -max [read len] -p <PAIRS file> -r <regions.bed> -w [window size]"
	echo ""
	echo "    -a            min read length"
	echo "    -b            max read length"
	echo "    -p              path to PAIRS.BED file"
	echo "    -r              regions of interest in BED format"
	echo "    -o              output path PREFIX"
	echo "    -w              size of desired window in bp. Single positive integer will result in symmetric padding"
	echo ""	
	exit 1
}

# if [ $# -ne 5 ]
# then
# 	display_usage
# 	exit 1
# fi

MIN=0
MAX=1000
PAIRSBED_FILE=""
REGIONS=""
OUTPREFIX=""
WIN=100

while getopts a:b:p:r:o:w:h flag; do
	case $flag in
	a) MIN=$OPTARG;;
	b) MAX=$OPTARG;;
	p) PAIRSBED_FILE=$OPTARG;;
	r) REGIONS=$OPTARG;;
	o) OUTPREFIX=$OPTARG;;
	w) WIN=$OPTARG;;
	h) display_usage;;
	esac
done

# Generate regions file with specified window parameters
echo "Determining windows around regions of interest"
WIN_REGIONS=$(basename "$REGIONS")
WIN_REGIONS=windowed.${WIN_REGIONS%.*}.tmp
${SCRIPTS_DIR}/bed_coord_window.pl -b $REGIONS -o $WIN_REGIONS -w $WIN


# Get PAIRS.BED file positions (fragments) that fall within regions of interest
echo "Getting reads that fall within windowed regions of interest"
REGION_READS=$(basename "$PAIRSBED_FILE")
REGION_READS=${REGION_READS}.tmp
${BEDOPS_DIR}/"bedops" --element-of -1 $PAIRSBED_FILE $WIN_REGIONS > $REGION_READS

# Process PAIRS.BED FILE and get oriented ends
${SCRIPTS_DIR}/average_plot_ends.pl -p $REGION_READS -r $WIN_REGIONS -o $OUTPREFIX -min $MIN -max $MAX -strand
