#!/bin/bash
. ./pairing.config

echo $2
TOTAL_MEM=$(( MEMORY_LAST*CORES_LAST/1000 ))

seqtk seq -a "$1" | lastal -m 100 -v -f BlastTab -P0 -i"$TOTAL_MEM"G -u 2 -Q 0 $PR2_LAST_DB | awk '$0 !~ /^#/ && $3 >= 95 && $4 >= 85' > mapped_cells_18s/"$2".last