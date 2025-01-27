#!/bin/bash

FILE_LIST=$1      # The subset file (barcode_subset_*)
THREADS=$2        # Number of threads
ASSEMBLY_DIR=$3   # Base directory for assembly

# Get the line corresponding to the current array task ID
line=$(sed -n "${LSB_JOBINDEX}p" $FILE_LIST)

# Extract the barcode and FASTQ file from the line
barcode=$(echo "$line" | cut -f1)
FASTQ=$(echo "$line" | cut -f2)

# Process the barcode and FASTQ
BARCODE_NO_WHITESPACE=$(echo "$barcode" | tr -d '[:space:]')
FASTQ_BASENAME=$(basename "$FASTQ")
BARCODE_FULL="$BARCODE_NO_WHITESPACE.$FASTQ_BASENAME"

if [ ! -f "$BARCODE_FULL/assembly/transcripts.fasta" ]; then
    mkdir -p "$BARCODE_FULL/assembly"

    for fq1 in $BARCODE_FULL/*L001*.filtered_trimmed_trimmed.fq.gz
    do
        fq2=$(echo ${fq1/L001/L002}) # Match the second pair
        echo $fq1
        if [ -f "$fq2" ]; then
            echo "Running SPAdes for $BARCODE_NO_WHITESPACE..."
            rnaspades.py --s1 "$fq1" --s2 "$fq2" -k 21,33,55 -t $THREADS -o "$BARCODE_FULL/assembly"
            echo "Assembly completed for $BARCODE_NO_WHITESPACE!"
        else
            echo "Warning: Paired file $fq2 not found for $fq1."
        fi
    done
else
    echo "Assembly already exists for $BARCODE_NO_WHITESPACE. Skipping..."
fi

