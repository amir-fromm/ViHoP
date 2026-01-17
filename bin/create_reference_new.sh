#!/bin/sh

module load bcl2fastq2 CellRanger/8.0.0

. ./pairing.config

## 2024-11-03
REFNAME=$1
FASTA=$2
GTF=$3


TOTAL_MEM=$(( CDHIT_CORES * CDHIT_MEM / 1000 ))

if ! ls "$DB_VIRUS_OUTPUT.$DEDUP.fasta" 1> /dev/null 2>&1
then
    sleep 2
    job_id=$(bsub -q $QUEUE_SHORT -J $WD/log/create_referece. -o $WD/log/create_referece.%J.o -e $WD/log/create_referece.%J.e  -n $CDHIT_CORES -R "rusage[mem=$CDHIT_MEM] span[ptile=$CDHIT_CORES]" cellranger mkref --genome=$REFNAME --fasta=$FASTA --nthreads $CDHIT_CORES --genes=$GTF --memgb=$TOTAL_MEM | awk '{print $2}' | tr -d '<>')
    job_ids+=($job_id)  # Store the job ID

    # Wait for all jobs to finish
    for id in "${job_ids[@]}"
    do
        while bjobs -noheader $id 2>/dev/null | grep -q 'PEND\|RUN\|UNKWN\|SSUSP'; do
            sleep 10  # Wait for 10 seconds before checking again
        done # Done checking status of jobs 
    done # The jobs are done
fi
