#!/bin/sh

module load bcl2fastq2 CellRanger/8.0.0

## 2024-11-03
REFNAME=$1
FASTA=$2
GTF=$3

cellranger mkref --genome=$REFNAME --fasta=$FASTA --genes=$GTF --memgb=16
