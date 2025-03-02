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
    # Detect all unique L00X lane identifiers
    read_sets=($(find "$BARCODE_FULL" -type f -name "*_L00*_*.filtered_trimmed_trimmed.fq.gz" | sed -E 's/.*_L00([0-9])_.*/L00\1/' | sort -u))

    # Debug: Show detected read sets
    echo "Detected read sets: ${read_sets[@]}"

    # Initialize SPAdes command
    SPADES_CMD="rnaspades.py -k 21,33,55 -t $THREADS -o $BARCODE_FULL/assembly"

    # Track argument index
    index=1

    for read_set in "${read_sets[@]}"; do
        # Find all FASTQ files matching this read set
        fq_files=($(find "$BARCODE_FULL" -type f -name "*${read_set}_*.filtered_trimmed_trimmed.fq.gz"))

        for fq in "${fq_files[@]}"; do
            if [ -f "$fq" ]; then
                SPADES_CMD+=" --s$index $fq"
                ((index++))
            else
                echo "Warning: Missing file for read set $read_set!"
            fi
        done
    done

    echo "Running SPAdes for $BARCODE_NO_WHITESPACE..."
    echo "$SPADES_CMD"  # Debugging: Print final SPAdes command
    eval "$SPADES_CMD"
    echo "Assembly completed for $BARCODE_NO_WHITESPACE!"

else
    echo "Assembly already exists for $BARCODE_NO_WHITESPACE. Skipping..."
fi

