#!/bin/bash
DIR=$1
ARRAY_JOB=$2
MEM=$3

mkdir -p log

ls $DIR > list_files.txt

split --numeric-suffixes -l 1000 list_files.txt list_files_subset

job_ids=()
for file_list in list_files_subset*
do
    number=$(bjobs -w | wc -l) # get the numbers of jobs running
    while [ "$number" -gt 400 ];do 
       echo "$number"
       sleep 1
       number=$(bjobs -w | wc -l)
    done
    total_files=$(wc -l < $file_list)
    job_id=$(bsub -n 1 -R "rusage[mem=$MEM]" -J "my_array_job[1-$total_files]" -o "log/%J_%I.o" -e "log/%J_%I.e" "$ARRAY_JOB $file_list" | awk '{print $2}' | tr -d '<>')
    job_ids+=($job_id)  # Store the job ID
done


# Wait for all jobs to finish
for id in "${job_ids[@]}"
do
    sleep 0.5
    while bjobs -noheader $id 2>/dev/null | grep -q 'PEND\|RUN\|UNKWN\|SSUSP'; do
        echo "sleeping"
        sleep 10  # Wait for 10 seconds before checking again
    done # Done checking status of jobs 
done # The jobs are done
job_ids=()
