#!/bin/bash

# Define the folder path
FOLDER_PATH=$1  

# Define the file patterns to delete
PATTERNS=(
    "*.filtered.fastq.gz"
    "*.filtered_trimmed.fq.gz"
    "*_trimming_report.txt"
    "*_trimmed_fastqc.html"
    "*_trimmed_fastqc.zip"
)

# Iterate through the patterns and delete matching files

for pattern in "${PATTERNS[@]}"; do
    find "$FOLDER_PATH" -type f -name "$pattern" -exec rm -v {} +
done


