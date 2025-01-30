#!/bin/sh
. ../../pairing.config

## Created 2025-01-29
## Extract reads from a fastq file based on their WD

echo "Pulling sequence headers by cell barcode"

# Input arguments
BARCODE_FILE=$1 # barcodes_and_fastq_paths.txt
MEMORY=$MEM_EXTRACT  # Initial memory allocation per job
MAX_JOBS=MAX_JOBS=$(( MAX_MEM /  MEM_EXTRACT ))
LOG_DIR="$WD/log_extract"

# Retry logic with a maximum iteration limit
max_retries=10  # Set a maximum number of retries
retry_count=0    # Initialize retry counter

mkdir -p $LOG_DIR

# Split barcode file for array jobs (500 jobs per batch)
split --numeric-suffixes -l 500 "$BARCODE_FILE" barcode_subset_

# Submit array jobs
job_ids=()
for subset_file in barcode_subset_*; do
    total_barcodes=$(wc -l < "$subset_file")
    echo "Submitting jobs for $subset_file with $total_barcodes tasks"

   job_id=$(bsub -n 1 -R "rusage[mem=$MEMORY]" \
    -J "extract_reads[1-$total_barcodes]" \
    -o "$LOG_DIR/%J_%I.o" -e "$LOG_DIR/%J_%I.e" \
    "sh $WD/bin/extract_reads_cell.sh \$(sed -n \${LSB_JOBINDEX}p $subset_file)" \
    | awk '{print $2}' | tr -d '<>')
done

# Wait for all jobs to finish
for id in "${job_ids[@]}"
do
    while bjobs -noheader $id 2>/dev/null | grep -q 'PEND\|RUN'; do
        sleep 10  # Wait for 10 seconds before checking again
    done # Done checking status of jobs 
done # The jobs are done


# Monitor and retry failed jobs if necessary
while :; do
    echo "Monitoring jobs..."
    sleep 10  # Check every 10 seconds
    failed_jobs=$(grep -l "completely empty" "$LOG_DIR"/*.e)
    if [ -z "$failed_jobs" ]; then
        echo "All jobs completed successfully!"
        break
    else
        retry_count=$((retry_count + 1))  # Increment retry counter

        if [ "$retry_count" -gt "$max_retries" ]; then
            echo "Maximum retry limit ($max_retries) reached. Exiting."
            exit 1
        fi
        echo "Some jobs failed. Retrying..."
        MEMORY=$((MEMORY * 2))  # Double memory for retries
        for failed_job_log in $failed_jobs; do
            job_index=$(basename "$failed_job_log" | cut -d'_' -f2 | cut -d'.' -f1)
            barcode=$(sed -n "${job_index}p" "$BARCODE_FILE")
            echo "Retrying job for barcode: $barcode"
            bsub -n 1 -R "rusage[mem=$MEMORY]" \
                -J "extract_retry" \
                -o "$LOG_DIR/%J.o" -e "$LOG_DIR/%J.e" \
                "sh $WD/bin/extract_reads_cell.sh $barcode"
        done
    fi
done