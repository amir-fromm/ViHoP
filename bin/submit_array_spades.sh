#!/bin/bash

. ./pairing.config

BARCODES=$1         # Tab-delimited file with barcodes and FASTQ paths
ARRAY_JOB=$2        # Script to execute per array job
MEM=$MEM_SPADES              # Initial memory allocation per job
THREADS=$(( CORES_SPADES * 2))           # Number of threads for SPAdes
ASSEMBLY_DIR="assemble_first_cells/all_cells"
BATCH_SIZE=500      # Number of jobs per batch
MAX_RUNNING_JOBS=$(( MAX_MEM / ( MEM_SPADES*CORES_SPADES ) )) # Maximum number of concurrent jobs
MAX_ITER=5  # Set the maximum number of iterations


cd assemble_first_cells/all_cells

mkdir -p $WD/log_spades $WD/tmp

# Function to submit jobs in batches
submit_jobs() {
    local barcode_file=$1
    local memory=$2
    local job_prefix=$3

    # Split the barcode file into smaller subsets for batching
    split --numeric-suffixes -l $BATCH_SIZE $barcode_file barcode_subset_

    job_ids=()
    for file_list in barcode_subset_*; do
        # Wait until the number of running jobs is below the maximum allowed
        number=$(bjobs -w | wc -l)
        while [ "$number" -gt $MAX_RUNNING_JOBS ]; do
            echo "Current running jobs: $number. Waiting..."
            sleep 10
            number=$(bjobs -w | wc -l)
        done

        # Submit the array job
        total_files=$(wc -l < $file_list)
        job_id=$(bsub -n $CORES_SPADES -R "rusage[mem=$memory]" \
            -J "${job_prefix}[1-$total_files]" \
            -o "$WD/log_spades/%J_%I.o" -e "$WD/log_spades/%J_%I.e" \
            "$ARRAY_JOB $file_list $THREADS $ASSEMBLY_DIR" | awk '{print $2}' | tr -d '<>')

        job_ids+=($job_id)
    done

    # Wait for all jobs to finish
    for id in "${job_ids[@]}"; do
        sleep 0.5
        while bjobs -noheader $id 2>/dev/null | grep -q 'PEND\|RUN\|UNKWN\|SSUSP'; do
            echo "Waiting for job $id to finish..."
            sleep 10
        done
    done

    # Cleanup temporary subsets
    rm barcode_subset_*
}

# Initial job submission
submit_jobs $BARCODES $MEM "spades_array_job"

# Retry loop for failed jobs
retry_needed=false
ITER=1

while [ $ITER -le $MAX_ITER ]; do
    echo "Checking for failed jobs (iteration $ITER)..."

    failed_jobs_file="$WD/tmp/failed_jobs_iteration_${ITER}.txt"
    > "$failed_jobs_file"  # Clear the file before use

    # **Check logs for memory failures ONLY**
    grep -l "reaching LSF memory usage limit" "$WD/log_spades/"* > "$WD/tmp/memory_limit_logs.txt" || true

    if [ -s "$WD/tmp/memory_limit_logs.txt" ]; then
        echo "Jobs failed due to memory issues. Collecting affected job indices..."
        
        # Extract job indices from log filenames
        grep -ohP "(?<=spades_array_job\[)\d+" "$WD/tmp/memory_limit_logs.txt" | sort -u > "$WD/tmp/memory_failed_indices.txt"

        # Find the corresponding barcodes for these failed indices
        while IFS= read -r index; do
            sed -n "${index}p" "$BARCODES" >> "$failed_jobs_file"
        done < "$WD/tmp/memory_failed_indices.txt"
    fi


    # **Retry failed jobs if there are any**
    if [ -s "$failed_jobs_file" ]; then
        retry_needed=true
        ((ITER++))
        rm $WD/log_spades/*
        if [ $ITER -gt $MAX_ITER ]; then
            echo "Maximum retries reached ($MAX_ITER). Exiting..."
            exit 1
        fi

        echo "Retrying failed jobs with increased memory (iteration $ITER)..."
        MEM=$((MEM * 2))  # Double the memory for retries
        submit_jobs "$failed_jobs_file" "$MEM" "spades_retry_job_$ITER"
        
        # Remove old log files only for failed jobs to avoid unnecessary re-runs
        xargs rm -f < "$WD/tmp/memory_limit_logs.txt"
    else
        echo "All jobs completed successfully!"
        exit 0
    fi
done

if [ "$retry_needed" = true ]; then
    echo "Some jobs failed after $MAX_ITER iterations. Check logs for more details."
fi


# Cleanup
rm -f $WD/tmp/memory_limit_logs.txt $WD/tmp/memory_failed_indices.txt $WD/tmp/failed_jobs_*.txt
