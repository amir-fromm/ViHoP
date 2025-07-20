
#!/bin/sh
. ./pairing.config
mkdir -p $WD/log
mkdir -p $WD/mapped_cells_virus

module load seqtk

LAST_DB=$1

job_ids=()

# Define the output file inside the job

out_file="$WD/mapped_cells_virus/last_candidates_all.last"


# Merge results of lastal mepping
(cat $WD/mapped_cells_virus/last_candidates_all.last $WD/mapped_cells_virus/last_viruses.last | grep -v '^#') > $WD/mapped_cells_virus/last_combined.last

# Summarize results (check if best match is to viruses)
LEN=$(cat $WD/mapped_cells_virus/last_combined.last | wc -l)
$WD/bin/summarize_lastal_one_large.py --chunks 100000 --size_table $LEN --output_path $WD/mapped_cells_virus/edited_last_combined.tsv --last_path $WD/mapped_cells_virus/last_combined.last

awk -F'\t' '!/^#/ && $14 == "Viruses" {OFS="\t"; print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}'  $WD/mapped_cells_virus/edited_last_combined.tsv >  $WD/mapped_cells_virus/summary_last_combined.tsv
