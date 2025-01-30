#!/bin/sh

if [ "$1" == "-h" ]; then
    echo "Usage: $0 <barcode> <FASTQ file>"
    exit 0
fi

barcode=$1
FASTQ=$2
FASTQB=$(basename "$FASTQ")

BARCODE_NO_WHITESPACE="$(echo -e "${barcode}" | tr -d '[:space:]')"
BARCODE_FULL=$BARCODE_NO_WHITESPACE.$FASTQB
mkdir -p "$BARCODE_FULL"
cd "$BARCODE_FULL"

if [ ! -f *.filtered_trimmed_trimmed.fq.gz ]; then
    # Create a single headers.temp file for all lanes to avoid redundant operations
    headers_file="$BARCODE_FULL.headers.temp"
    for FASTQR1 in "$FASTQ"/*_R1_*.fastq.gz; do
        echo "Processing $FASTQR1"
        zgrep -B1 "$BARCODE_NO_WHITESPACE" "$FASTQR1" | grep '@' | cut -d' ' -f1 | tr -d '@' >> "$headers_file"
    done

    # Remove duplicates from headers.temp to speed up seqtk subseq
    sort -u "$headers_file" -o "$headers_file"

    # Process each FASTQR2 (R2 files) based on the headers
    for FASTQR1 in "$FASTQ"/*_R1_*.fastq.gz; do
        FASTQR2=${FASTQR1/R1/R2}
        echo "Filtering sequences from $FASTQR2"
        f=$(basename "$FASTQR2" .fastq.gz)
        OUT="${f}_${BARCODE_FULL}.filtered.fastq.gz"
        seqtk subseq "$FASTQR2" "$headers_file" | gzip > "$OUT"
    done

    # Remove headers.temp file after use
    rm "$headers_file"

    # Trimming step 1
    echo "Trimming sequences (first pass)..."
    for fastq in *".filtered.fastq.gz"; do
        if [ -f "$fastq" ]; then
            trim_galore --phred33 -j 8 --length 36 -q 5 --stringency 1 --gzip --fastqc -e 0.1 "$fastq"
            sleep 1
        fi
    done

    # Trimming step 2 (polyA trimming)
    echo "Trimming sequences (polyA removal)..."
    for fastq in *"trim"*".fq.gz"; do
        if [ -f "$fastq" ]; then
            trim_galore --polyA -j 1 --length 36 --gzip "$fastq"
            sleep 1
        fi
    done
fi
