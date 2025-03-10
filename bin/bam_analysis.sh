#!/bin/bash

#module load spider
module load SAMtools
set -e

BAM_FILE=$1
OUT_FOLDER=$2
SAMPLE=$3
THREADS=$4

mkdir -p $OUT_FOLDER

MAPPED_BAM=${BAM_FILE%.bam}.mapped.bam

# get total reads per barcode for all barcodes (with at least one read)
#echo "Getting total reads per barcode"
#samtools view --threads $THREADS $BAM_FILE | grep -P 'CB:Z:' | awk '{for(i=1;i<=NF;i++) if($i ~ /^CB:Z:/) {split($i,a,":"); print a[3]}}' | sort --parallel=$THREADS | uniq -c > $OUT_FOLDER/$SAMPLE.reads_per_barcode.txt

# extract only mapped reads subset from bam file
echo "Extract only mapped reads from bam file"
samtools view --threads $THREADS -b -F 4 $BAM_FILE  > $MAPPED_BAM
samtools index $MAPPED_BAM

# get UMI counts per barcode
echo "Getting UMI counts per barcode"
samtools view --threads $THREADS $MAPPED_BAM | awk '/xf:i:25/ && /CB:Z:/ {for(i=1;i<=NF;i++) if($i ~ /^CB:Z:/) {split($i,a,":"); print a[3]}}' | sort --parallel=$THREADS | uniq -c > $OUT_FOLDER/$SAMPLE.mapped_UMIs_per_barcode.txt

# get mapped reads counts
echo "Getting mapped reads counts"
samtools view --threads $THREADS $MAPPED_BAM | awk '($5 == 255 || $1 ~ /^@/) && /CB:Z:/ {for(i=1;i<=NF;i++) if($i ~ /^CB:Z:/) {split($i,a,":"); print a[3]}}' | sort --parallel=$THREADS | uniq -c > $OUT_FOLDER/$SAMPLE.mapped_subset.txt

if [[ -s $MAPPED_BAM  ]]
then
    rm $BAM_FILE
fi
echo "Done!"