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
    job_ids+=("$job_id")  # <-- ADD JOB ID TO ARRAY

done

# Wait for all jobs to finish
echo "Waiting for all jobs to finish..."
while bjobs -J "extract_reads" 2>/dev/null | grep -q 'PEND\|RUN\|UNKWN\|SSUSP'; do
    sleep 10
done
echo "All jobs completed!"


# Monitor and retry failed jobs if necessary

while :; do
    echo "Monitoring jobs..."
    sleep 10  # Check every 10 seconds

    # Identify jobs that failed due to memory limit
    failed_jobs=$(grep -l "reaching LSF memory usage limit" "$LOG_DIR"/*.e || true)

    if [ -z "$failed_jobs" ]; then
        echo "All jobs completed successfully!"
        break
    fi

    retry_count=$((retry_count + 1))  # Increment retry counter

    if [ "$retry_count" -gt "$max_retries" ]; then
        echo "Maximum retry limit ($max_retries) reached. Exiting."
        exit 1
    fi

    num_failed_jobs=$(echo "$failed_jobs" | wc -l)
    if [ "$num_failed_jobs" -eq 0 ]; then
        echo "No jobs failed due to memory. Exiting."
        break
    fi

    # Calculate new memory per job
    new_memory=$(($MAX_MEM / num_failed_jobs))
    
    # Ensure we do not exceed total memory
    if [ "$new_memory" -le 0 ]; then
        echo "Not enough memory available to retry jobs. Exiting."
        exit 1
    fi

    echo "Some jobs failed due to memory. Retrying with $new_memory MB per job..."

    # Retry only failed jobs
    for failed_job_log in $failed_jobs; do
        job_index=$(basename "$failed_job_log" | cut -d'_' -f2 | cut -d'.' -f1)
        barcode=$(sed -n "${job_index}p" "$BARCODE_FILE")

        echo "Retrying job for barcode: $barcode with $new_memory MB"

        bsub -n 1 -R "rusage[mem=$new_memory]" \
            -J "extract_retry_$retry_count" \
            -o "$LOG_DIR/%J.o" -e "$LOG_DIR/%J.e" \
            "sh $WD/bin/extract_reads_cell.sh $barcode"
    done
done