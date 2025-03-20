#!/bin/bash
. ./pairing.config

echo $2
TOTAL_MEM=$(( MEMORY_LAST*CORES_LAST/1000 ))

seqtk seq -a "$1" | lastal -m 100 -v -f BlastTab -P0 -i"$TOTAL_MEM"G -u 2 -Q 0 $PR2_LAST_DB | awk -v id="$LAST_IDENTITY" -v len="$LAST_LENGTH" '$0 !~ /^#/ && $3 >= id && $4 >= len' > mapped_cells_18s/"$2".last