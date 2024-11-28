#!/bin/sh


# Get the array index from LSB_JOBINDEX
index="$LSB_JOBINDEX"

file_list=$1
# Read the file names from file_list.txt into an array
mapfile -t files < $file_list

# Get the file name using the array index
input_file="${files[$index - 1]}"  # Subtract 1 to match 0-based indexing
CELL=$(basename $input_file)
echo $CELL
sed "s/>.*/&_$CELL/" $CELL/assembly/transcripts.fasta > $CELL/assembly/transcripts.edit.fasta

