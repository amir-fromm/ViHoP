
#!/bin/sh
. ./pairing.config
mkdir -p $WD/log
mkdir -p $WD/mapped_cells_virus

module load seqtk

job_ids=()

# Define the output file inside the job
out_file="$WD/mapped_cells_virus/last_viruses.last"

if [ ! -f "$out_file" ]
then
    echo "Running virus lastal"
    sleep 1
    # Submit the job and capture job ID
    job_id=$(bsub -q short \
    -J lastal_virus \
    -o "$WD/log/lastal_virus.%J.o" \
    -e "$WD/log/lastal_virus.%J.e" \
    -n $CORES_SEARCH \
    -R "rusage[mem=$MEMORY_SEARCH] span[ptile=$CORES_SEARCH]" \
    "lastal -m 100 -f BlastTab -P 0 -i14G -u 2 -v -Q 0 -F 15 $LAST_DB $WD/assemble_first_cells/all_cells.transcripts.final.fasta > $out_file" \
    | awk '{print $2}' | tr -d '<>')

    job_ids+=("$job_id")  # Store the job ID

    # Wait for all jobs to finish
    for id in "${job_ids[@]}"
    do
        while bjobs -noheader "$id" 2>/dev/null | grep -q 'PEND\|RUN\|UNKWN\|SSUSP'; do
            sleep 10
        done
    done
fi
# extract only matching contigs to verify them agaisnt the general database
echo "Extract matching contigs to verify them agaisnt the general database"
awk '!/^#/ {print $1}' $WD/mapped_cells_virus/last_viruses.last | sort -u > $WD/mapped_cells_virus/last_viruses_candidates.txt

seqtk subseq $WD/assemble_first_cells/all_cells.transcripts.final.fasta $WD/mapped_cells_virus/last_viruses_candidates.txt > $WD/mapped_cells_virus/last_viruses_candidates.fasta

# lastal of the matching single-cell assembled contigs against a general protein database
job_ids=()

out_file="$WD/mapped_cells_virus/last_candidates_all.last"
if [ ! -f "$out_file" ]
then
    echo "Running general lastal"
    sleep 1
    job_id=$(bsub -q short \
    -J $WD/log/lastal_refseq. \
    -o $WD/log/lastal_refseq.%J.o \
    -e $WD/log/lastal_refseq.%J.e  \
    -n $CORES_SEARCH \
    -R "rusage[mem=$MEMORY_SEARCH] span[ptile=$CORES_SEARCH]" \
    "lastal -m 100 -f BlastTab -P 0 -i32G -u 2 -v -Q 0 -F 15 /home/labs/vardi/Collaboration/databases/Aylward_DB_2024/refSeq223_Virus_appended_masked_lin $WD/mapped_cells_virus/last_viruses_candidates.fasta > $out_file" \
    | awk '{print $2}' | tr -d '<>')
    job_ids+=($job_id)  # Store the job ID

    # Wait for all jobs to finish
    for id in "${job_ids[@]}"
    do
        while bjobs -noheader $id 2>/dev/null | grep -q 'PEND\|RUN\|UNKWN\|SSUSP'; do
            sleep 10  # Wait for 10 seconds before checking again
        done # Done checking status of jobs 
    done # The jobs are done
fi
##

# Merge results of lastal mepping
#echo "Merge results of lastal mepping"
#(cat $WD/mapped_cells_virus/last_candidates_all.last $WD/mapped_cells_virus/last_viruses.last | grep -v '^#') > $WD/mapped_cells_virus/last_combined.last

# Summarize results (check if best match is to viruses)
LEN=$(cat $WD/mapped_cells_virus/last_candidates_all.last | wc -l)
$WD/bin/summarize_lastal_one_large.py --chunks 100000 --size_table $LEN --output_path $WD/mapped_cells_virus/edited_last_all.tsv --last_path $WD/mapped_cells_virus/last_candidates_all.last
echo "Print summarized table"
awk -F'\t' '
  !/^#/ && ($14 == "Bacteria" || $14 == "Eukaryotes" || $14 == "Others" || $14 == "Phages") {
    if (!seen[$1]++) print $1
  }
' $WD/mapped_cells_virus/edited_last_all.tsv > $WD/mapped_cells_virus/false_positives.list
#awk -F'\t' '!/^#/ && $14 == "Viruses" {OFS="\t"; print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}'  $WD/mapped_cells_virus/edited_last_combined.tsv >  $WD/mapped_cells_virus/summary_last_combined.tsv
(cat $WD/mapped_cells_virus/last_candidates_all.last $WD/mapped_cells_virus/last_viruses.last | grep -v '^#') > $WD/mapped_cells_virus/last_combined.last

# Summarize results (check if best match is to viruses)
LEN=$(cat $WD/mapped_cells_virus/last_combined.last | wc -l)
$WD/bin/summarize_lastal_one_large.py --chunks 100000 --size_table $LEN --false_positives $WD/mapped_cells_virus/false_positives.list --output_path $WD/mapped_cells_virus/edited_last_combined.tsv --last_path $WD/mapped_cells_virus/last_combined.last

FINAL_FILE=$WD/mapped_cells_virus/summary_last_combined.tsv
if [ ! -f "$FINAL_FILE" ]
    echo "Virus last fail!"
awk -F'\t' '!/^#/ && $14 == "Viruses" {OFS="\t"; print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}'  $WD/mapped_cells_virus/edited_last_combined.tsv >  $WD/mapped_cells_virus/summary_last_combined.tsv
