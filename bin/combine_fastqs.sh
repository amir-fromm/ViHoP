#!/bin/bash

# Usage: combine_fastqs.sh folder_list.txt

folder_list="$1"

# Get the folder corresponding to this array index
folder=$(sed -n "${LSB_JOBINDEX}p" "$folder_list")

if [[ -z "$folder" ]]; then
    echo "No folder found for index $LSB_JOBINDEX"
    exit 1
fi

# Remove trailing slash if any
folder="${folder%/}"
CELL=$(basename "$folder")

echo "Processing folder: $folder"

# Combine all matching fastq.gz files
cat "$folder"/*.filtered_trimmed_trimmed.fq.gz > "$folder/$CELL.combined.fq.gz"
