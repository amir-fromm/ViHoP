#!/bin/bash

. ./pairing.config

BARCODES=$1         # Tab-delimited file with barcodes and FASTQ paths
ARRAY_JOB=$2        # Script to execute per array job
MEM=$MEM_SPADES              # Initial memory allocation per job
THREADS=$(( CORES_SPADES * 2))           # Number of threads for SPAdes
ASSEMBLY_DIR="assemble_first_cells/all_cells"
BATCH_SIZE=500      # Number of jobs per batch
MAX_RUNNING_JOBS=$(( MAX_MEM / ( MEM_SPADES*CORES_SPADES ) )) # Maximum number of concurrent jobs
ITER=1              # Retry iteration count


cd assemble_first_cells/all_cells

mkdir -p log_array tmp

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
            -o "log_array/%J_%I.o" -e "log_array/%J_%I.e" \
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
retry_needed=true
while $retry_needed && [ "$ITER" -lt 10 ]; do
    retry_needed=false
    failed_jobs_file="tmp/failed_jobs_$ITER.txt"
    > $failed_jobs_file  # Empty the file

    echo "Checking for failed jobs (iteration $ITER)..."

    # Check logs for memory errors
    grep -l "reaching LSF memory usage limit" log_array/* > tmp/memory_limit_logs.txt || true
    if [ -s tmp/memory_limit_logs.txt ]; then

        echo "Jobs failed due to memory issues. Collecting affected jobs..."
        grep -ohP "(?<=spades_array_job\[)\d+" tmp/memory_limit_logs.txt | sort -u > tmp/memory_failed_indices.txt

        # Find the corresponding barcodes for these indices
        while IFS= read -r index; do
            sed -n "${index}p" $BARCODES >> $failed_jobs_file
        done < tmp/memory_failed_indices.txt
    fi

    # Fallback: Check for missing output files
    while IFS=$'\t' read -r barcode FASTQ; do
        barcode_no_whitespace=$(echo "$barcode" | tr -d '[:space:]')
        fastq_basename=$(basename "$FASTQ")
        barcode_full="$ASSEMBLY_DIR/$barcode_no_whitespace.$fastq_basename"

        if [ ! -f "$barcode_full/assembly/transcripts.fasta" ]; then
            echo -e "$barcode\t$FASTQ" >> $failed_jobs_file
        fi
    done < $BARCODES

    # Retry failed jobs if there are any
    if [ -s $failed_jobs_file ]; then
        retry_needed=true
        ((ITER++))
        echo "Retrying failed jobs with increased memory (iteration $ITER)..."
        MEM=$((MEM * 2))  # Double the memory for retries
        submit_jobs $failed_jobs_file $MEM "spades_retry_job_$ITER"
        rm -f log_array/*
    fi
done

if ! $retry_needed; then
    echo "All jobs completed successfully!"
else
    echo "Some jobs failed after $ITER iterations. Check logs for more details."
fi

# Cleanup
rm -f tmp/memory_limit_logs.txt tmp/memory_failed_indices.txt tmp/failed_jobs_*.txt
