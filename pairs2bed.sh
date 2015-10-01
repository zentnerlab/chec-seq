#!/bin/bash

#Siva Kasinathan (Henikoff Lab, FHCRC): skasin@uw.edu

#Converts PAIRS file into UCSC-compatible BED
#Dependencies: BEDOPS


#PAIRS format:
# Column	Value
# 0			Chromosome
# 1			Strand
# 2 		Start
# 3			End
# 4			Map quality
# 5 		Length

#PAIRS.BED format:
# Column	Value
# 0			Chromosome
# 1			Start
# 2			End
# 3			Length
# 4			Map quality
# 5			Strand

BEDOPS_DIR="/home/skasinat/bedops/bin/"

display_usage(){
	echo ""
	echo "Usage: ./pairs2bed.sh <PAIRS file>"
	echo ""
	echo "    <PAIRS file> path to PAIRS file"
	echo ""
}

if [ $# -ne 1 ]
then
	display_usage
	exit 1
fi

PAIRS=$1
OUT_FILE=$(basename "$PAIRS")
OUT_FILE="${OUT_FILE%.*}"".pairs.bed"
OUTPUT=$PWD"/"$OUT_FILE

echo ""
echo "Converting .pairs file to .pairs.bed file"
awk '{print $1 "\t" $3-1 "\t" $4 "\t" $6 "\t" $5 "\t" $2}' $PAIRS > $OUTPUT
# Note that PAIRS files appear to be indexed starting at 1 (i.e., first position
# in a chromosome has index 1). In contrast, UCSC BED files are 0 indexed.
# Therefore, 1 is subtracted from PAIRS start entries to effect conversion.

echo "Sorting .pairs.bed file"
$BEDOPS_DIR"sort-bed" $OUTPUT > temp
mv temp $OUTPUT
echo "Done"
echo ""
