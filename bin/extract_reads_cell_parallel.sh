#!/bin/sh
. ../../pairing.config

## Created 2022-05-02

# updated to other versions of modules in 22-8-24
# updated to a streamline protocol in 3-11-24

BARCODES=$1 # List of cell barcodes. E.g. cells cells with high probability of infection
job_ids=()
echo "Pulling sequence headers by cell barcode"

mkdir -p $WD/log_extract
MAX_JOBS=$(( MAX_MEM /  MEM_EXTRACT ))

while IFS=$'\t' read -r barcode FASTQ || [ -n "$barcode" ]
    do

        FASTQB=$(basename  $FASTQ)
        BARCODE_NO_WHITESPACE="$(echo -e "${barcode}" | tr -d '[:space:]')"
        BARCODE_FULL=$BARCODE_NO_WHITESPACE.$FASTQB

        if ! ls "$WD/assemble_first_cells/all_cells/$BARCODE_FULL/"*.filtered_trimmed_trimmed.fq.gz 1> /dev/null 2>&1; then

            echo "$BARCODE_FULL"
            number=$(bjobs -w | wc -l) # get the numbers of jobs running
            while [ "$number" -gt $MAX_JOBS ];do 
                echo "$number"
                sleep 1
                number=$(bjobs -w | wc -l)
            done
            
            job_id=$(bsub -q $QUEUE_SHORT -J $WD/log_extract/extract_reads_$barcode -o $WD/log_extract/extract_reads_$barcode.%J.o -e $WD/log_extract/extract_reads_$barcode.%J.e -n 1 -R "rusage[mem=$MEM_EXTRACT]" sh $WD/bin/extract_reads_cell.sh $barcode $FASTQ  | awk '{print $2}' | tr -d '<>')
            job_ids+=($job_id)  # Store the job ID
            
            sleep 0.01
        fi
    done  < $BARCODES

# Wait for all jobs to finish
for id in "${job_ids[@]}"
do
    while bjobs -noheader $id 2>/dev/null | grep -q 'PEND\|RUN'; do
        sleep 10  # Wait for 10 seconds before checking again
    done # Done checking status of jobs 
done # The jobs are done


