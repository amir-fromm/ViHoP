#!/bin/sh

. ../../pairing.config


barcode=$1
FASTQ=$2
CORES=$3

THREADS=$(( 2*CORES ))

FASTQB=$(basename  $FASTQ)

BARCODE_NO_WHITESPACE="$(echo -e "${barcode}" | tr -d '[:space:]')"
BARCODE_FULL=$BARCODE_NO_WHITESPACE.$FASTQB
echo "$BARCODE_FULL"
cd $BARCODE_FULL

mkdir -p assembly

for FQ in *L001*.filtered_trimmed_trimmed.fq.gz
do
    NAME="$(basename $FQ .fq.gz)"
    fq1=$FQ
    fq2="$(echo ${FQ/L001/L002})" 
    echo "Running spades"
    ## run script
    rnaspades.py --s 1 $fq1 --s 2 $fq2 -k 21 33 55 -t $THREADS -o assembly
    echo "Done!"
    sleep 1
done
    # rm *.fq.gz
    # rm *.filtered*

cd ..
