
#!/bin/bash
. ./pairing.config

input_list="$1"  # List of files
if [[ -z "$LSB_JOBINDEX" ]]; then
    echo "Error: LSB_JOBINDEX is not set."
    exit 1
fi

input_file=$(sed -n "${LSB_JOBINDEX}p" "$input_list")  # Get the file for this job index

if [[ ! -s "$input_file" ]]; then
    echo "Skipping empty or missing file: $input_file"
    exit 1
fi

filename="${input_file##*/}"
filename="${filename%.combined.fq.gz}"
TOTAL_MEM=$(( MEMORY_LAST*CORES_LAST/1000 ))

echo "Processing: $input_file"

seqtk seq -a "$input_file" | \
lastal -m 100 -v -f BlastTab -P0 -i"${TOTAL_MEM}G" -u 2 -Q 0 $PR2_LAST_DB | \
awk -v id="$LAST_IDENTITY" -v len="$LAST_LENGTH" '$0 !~ /^#/ && $3 >= id && $4 >= len' > \
mapped_cells_18s/"$filename".last
