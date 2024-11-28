#!/bin/bash

. ./pairing.config

MAX_JOBS=$(( MAX_MEM / ( MEMORY_LAST*CORES_LAST ) ))
job_ids=()
while IFS= read -r file || [[ -n "$file" ]]
do
    if [[ -s "$file" ]]
    then
        input_file=$file
        filename="${input_file##*/}"
        filename="${filename%.combined.fq.gz}"
        if ! ls "$WD/mapped_cells_18s/$filename"*.last 1> /dev/null 2>&1; then
            echo "Processing file: $filename"
            number=$(bjobs -w | wc -l) # get the numbers of jobs running
            while [ "$number" -gt $MAX_JOBS ];do 
                echo "$number"
                sleep 1
                number=$(bjobs -w | wc -l)
            done
            
            job_id=$( bsub -q $QUEUE_SHORT -J log_last/18s_lastal.$filename. -o log_last/18s_lastal.$filename.%J.o -e log_last/18s_lastal.$filename.%J.e -R "rusage[mem=$MEMORY_LAST] affinity[core(1):distribute=balance] span[ptile=$CORES_LAST]" -n $CORES_LAST bin/lastal_loop_18s_run.sh $input_file $filename | awk '{print $2}' | tr -d '<>')
            job_ids+=($job_id)
            sleep 0.01
        fi
    fi
    
done < "$1"
 

for id in "${job_ids[@]}"
do
    while bjobs -noheader $id 2>/dev/null | grep -q 'PEND\|RUN'; do
        sleep 10  # Wait for 10 seconds before checking again
    done # Done checking status of jobs 
done # The jobs are done




    