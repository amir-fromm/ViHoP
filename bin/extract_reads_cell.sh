#!/bin/sh

if [ "$1" == "-h" ]; then
    echo ""
    exit 0
fi

barcode=$1
FASTQ=$2
FASTQB=$(basename  $FASTQ)

BARCODE_NO_WHITESPACE="$(echo -e "${barcode}" | tr -d '[:space:]')"
BARCODE_FULL=$BARCODE_NO_WHITESPACE.$FASTQB
mkdir -p $BARCODE_FULL
cd $BARCODE_FULL

if [ ! -f *.filtered_trimmed_trimmed.fq.gz ]; then
    for FASTQR1 in $FASTQ/*_R1_*.fastq.gz # Iterate fastq files from all four lanes
        do
        echo $FASTQR1
        zcat $FASTQR1 | grep -B1 $BARCODE_NO_WHITESPACE | grep @ | tr -d @ | awk '{print $1;}' >> $FASTQB.headers.temp # Greb sequence headers by cell barcode recorded in R1
        FASTQR2=${FASTQR1/R1/R2} # Switch to R2
        echo "Filtering fastq sequences"
        echo $FASTQR2
        f="$(basename $FASTQR2 .fastq.gz)"
        OUT=$f"_"$BARCODE_FULL".filtered.fastq.gz"
        echo $OUT
        seqtk subseq $FASTQR2 $FASTQB.headers.temp | gzip > $OUT # Get sequences stored in R4 by headers
        rm $FASTQB.headers.temp
        done

        for fastq in *".fastq.gz"
            do
            if [ -f "$fastq" ];then
                echo "$fastq"
                trim_galore --phred33 -j 8 --length 36 -q 5 --stringency 1 --gzip --fastqc -e 0.1 $fastq
                sleep 1
            fi
            done

        for fastq in *"trim"*".fq.gz"
            do
            if [ -f "$fastq" ];then
                echo "$fastq"
                trim_galore --polyA -j 1 --length 36 --gzip $fastq
                sleep 1
            fi
            done
fi

