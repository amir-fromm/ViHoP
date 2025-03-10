#!/bin/sh

if [ "$1" == "-h" ]; then
    echo "Usage: $0 <barcode> <FASTQ file(s)>"
    exit 0
fi

barcode=$1
FASTQ_LIST=$2  # List of one or multiple FASTQ paths, comma-separated
FASTQ_PATHS=($(echo "$FASTQ_LIST" | tr ',' ' '))  # Convert comma-separated list into array

# Extract FASTQ base name from the first path
FASTQB=$(basename "${FASTQ_PATHS[0]}")

BARCODE_NO_WHITESPACE="$(echo -e "${barcode}" | tr -d '[:space:]')"
BARCODE_FULL="${BARCODE_NO_WHITESPACE}.${FASTQB}"  # Use full barcode format

mkdir -p "$BARCODE_FULL"
cd "$BARCODE_FULL" || exit 1

if [ ! -f *.filtered_trimmed_trimmed.fq.gz ]; then
    headers_file="$BARCODE_FULL.headers.temp"
    rm -f "$headers_file"  # Ensure no old headers file exists

    # Step 1: Collect sequence headers from all FASTQ paths
    for FASTQ in "${FASTQ_PATHS[@]}"; do
        for FASTQR1 in "$FASTQ"/*_R1_*.fastq.gz; do

            FASTQ_NAME=$(basename "${FASTQR1}")
            headers_file="${FASTQ_NAME}_headers.temp"  # Unique header file per R1 file
            rm -f "$headers_file"

            echo "Processing $FASTQR1"
            echo "Headers file: $headers_file"
            zgrep -B1 "$BARCODE_NO_WHITESPACE" "$FASTQR1" | grep '@' | cut -d' ' -f1 | tr -d '@' > "$headers_file"
            if [ ! -s "$headers_file" ]; then
                echo "Warning: No reads found for barcode $BARCODE_NO_WHITESPACE in $FASTQR1!"
                continue
            fi
            # Step 2: Find corresponding R2 file
            FASTQR2="${FASTQR1/_R1_/_R2_}"
            if [ -f "$FASTQR2" ]; then
                echo "Filtering sequences from $FASTQR2"
                f=$(basename "$FASTQR2" .fastq.gz)
                OUT="${f}_${BARCODE_FULL}.filtered.fastq.gz"

                # Extract reads using the correct header file
                seqtk subseq "$FASTQR2" "$headers_file" | gzip > "$OUT"
            else
                echo "Warning: Missing expected R2 file for $FASTQR1"
            fi

            rm -f "$headers_file"  # Clean up temp header file
        done
    done


    # Step 3: Trimming sequences
    for fastq in *".filtered.fastq.gz"; do
        if [ -f "$fastq" ]; then
            trim_galore --phred33 -j 8 --length 36 -q 5 --stringency 1 --gzip --fastqc -e 0.1 "$fastq"
            sleep 1
        fi
    done

    for fastq in *"trim"*".fq.gz"; do
        if [ -f "$fastq" ]; then
            trim_galore --polyA -j 1 --length 36 --gzip "$fastq"
            sleep 1
        fi
    done
fi
