#!/bin/sh

module load bcl2fastq2 CellRanger/8.0.0

if [ "$1" == "-h" ]; then
    echo "a script that runs through all fastq folders given a specific regex.
    first variable is reference
    second variable is regex for fastq path
    third variable is memory"
    exit 0
fi
mkdir -p log
REF=$1
FASTQDIR=$2
MEM=$3
CORES=$4
ID=$(basename $FASTQDIR)
DIRPATH="$(dirname "${REF}")"
echo "$DIRPATH"
### Edited: 2024-07-22 ###
for dir in $FASTQDIR*/; do
    echo $dir
    ID=$(basename $dir)
    bsub -q short -J log/cellranger_$ID -o log/cellranger_$ID%J.o -e log/cellranger_$ID%J.e  -n $CORES -R "rusage[mem=$MEM] span[ptile=$CORES]" cellranger count --id=$ID --transcriptome=$REF --fastqs=$dir --project=$ID --create-bam=true --expect-cells=10000 --chemistry=SC3Pv3  --nosecondary
done
