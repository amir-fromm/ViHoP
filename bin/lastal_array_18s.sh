#!/bin/bash

. ./pairing.config  # Load configuration variables

MAX_RETRIES=5  # Limit number of retries
MAX_BATCH_SIZE=2000  # Max job array size allowed
retry_count=1

input_list="$WD/tmp/combined_reads.txt"
failed_input_list="$WD/tmp/failed_combined_reads.txt"

# Function to submit jobs in batches
submit_jobs() {
    local list="$1"
    local memory="$2"
    local retry="$3"

    echo "Submitting job array for retry #$retry with memory: $memory MB"

    total_jobs=$(wc -l < "$list")
    batch_start=1

    while [ "$batch_start" -le "$total_jobs" ]; do
        batch_end=$(( batch_start + MAX_BATCH_SIZE - 1 ))
        if [ "$batch_end" -gt "$total_jobs" ]; then
            batch_end="$total_jobs"
        fi

        job_id=$(bsub -q $QUEUE_SHORT \
            -J "lastal_retry_${retry}[${batch_start}-${batch_end}]" \
            -o "log_last/18s_lastal_retry_${retry}.%I.o" \
            -e "log_last/18s_lastal_retry_${retry}.%I.e" \
            -R "rusage[mem=$memory] affinity[core(1):distribute=balance] span[ptile=$CORES_LAST]" \
            -n $CORES_LAST \
            bin/lastal_array_18s_run.sh "$list" | awk '{print $2}' | tr -d '<>')

        echo "Submitted batch: ${batch_start}-${batch_end}, Job ID: $job_id"

        batch_start=$(( batch_end + 1 ))
    done

    # Wait for jobs to complete
    echo "Waiting for jobs to complete..."
    while bjobs -J "lastal_retry_${retry}" 2>/dev/null | grep -q 'PEND\|RUN'; do
        sleep 30
    done
}

# Initial job submission
submit_jobs "$input_list" "$MEMORY_LAST" "$retry_count"

# Retry loop for memory failures
while [ $retry_count -lt $MAX_RETRIES ]; do
    echo "Checking for failed jobs..."
    
    # Capture failed job logs
    failed_jobs=($(grep -l "TERM_MEMLIMIT" log_last/*.o))
    rm -f log_last/*o log_last/*e

    if [ ${#failed_jobs[@]} -eq 0 ]; then
        echo "No failed jobs due to memory limits."
        break
    fi

    echo "Failed jobs detected: ${#failed_jobs[@]}"

    # Extract numeric job indices from logs
    failed_indices=()
    for job in "${failed_jobs[@]}"; do
        job_index=$(echo "$job" | sed -E 's/.*18s_lastal_retry_[0-9]+\.([0-9]+)\..*/\1/')
        if [[ "$job_index" =~ ^[0-9]+$ ]]; then
            echo "Extracted Job Index: $job_index"
            failed_indices+=("$job_index")
        else
            echo "Warning: Could not extract job index from $job"
        fi
    done

    if [ ${#failed_indices[@]} -eq 0 ]; then
        echo "No valid job indices found. Exiting."
        break
    fi

    # Create a new failed input list
    rm -f "$failed_input_list"
    for index in "${failed_indices[@]}"; do
        sed -n "${index}p" "$input_list" >> "$failed_input_list"
    done

    if [[ ! -s "$failed_input_list" ]]; then
        echo "Error: No failed input list created!"
        break
    fi

    # Double memory for retries
    MEMORY_LAST=$(( MEMORY_LAST * 2 ))
    retry_count=$(( retry_count + 1 ))

    echo "Retrying failed jobs (Attempt $retry_count)..."
    submit_jobs "$failed_input_list" "$MEMORY_LAST" "$retry_count"
done

echo "All retries completed."
