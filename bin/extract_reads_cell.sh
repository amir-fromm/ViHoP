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
            echo "Processing $FASTQR1"
            zgrep -B1 "$BARCODE_NO_WHITESPACE" "$FASTQR1" | grep '@' | cut -d' ' -f1 | tr -d '@' >> "$headers_file"
        done
    done

    sort -u "$headers_file" -o "$headers_file"  # Remove duplicates

    # Step 2: Extract matching reads from R2 files
    for FASTQ in "${FASTQ_PATHS[@]}"; do
        for FASTQR1 in "$FASTQ"/*_R1_*.fastq.gz; do
            FASTQR2="${FASTQR1/_R1_/_R2_}"  # Correct substitution for R2
            if [ -f "$FASTQR2" ]; then
                echo "Filtering sequences from $FASTQR2"
                f=$(basename "$FASTQR2" .fastq.gz)
                OUT="${f}_${BARCODE_FULL}.filtered.fastq.gz"  # Use full barcode format
                seqtk subseq "$FASTQR2" "$headers_file" | gzip > "$OUT"
            else
                echo "Warning: Missing expected R2 file for $FASTQR1"
            fi
        done
    done

    rm -f "$headers_file"  # Remove temp file after use

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
