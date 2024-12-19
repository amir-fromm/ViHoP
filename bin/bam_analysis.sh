
#module load spider
module load SAMtools

BAM_FILE=$1
OUT_FOLDER=$2
SAMPLE=$3
THREADS=$4

mkdir -p $OUT_FOLDER
# get total reads per barcode for all barcodes (with at least one read)
echo "Getting total reads per barcode"
#samtools view --threads $THREADS $BAM_FILE | grep CB:Z: | sed 's/.*CB:Z:\([ACGT]*\).*/\1/' | sort | uniq -c > $OUT_FOLDER/$SAMPLE.reads_per_barcode.txt
samtools view --threads $THREADS $BAM_FILE | grep -P 'CB:Z:' | awk '{for(i=1;i<=NF;i++) if($i ~ /^CB:Z:/) {split($i,a,":"); print a[3]}}' | sort --parallel=$THREADS | uniq -c > $OUT_FOLDER/$SAMPLE.reads_per_barcode.txt
# get UMI counts per barcode
echo "Getting UMI counts per barcode"
#samtools view --threads $THREADS $BAM_FILE | grep xf:i:25 | sed 's/.*CB:Z:\([ACGT]*\).*/\1/' | sort | uniq -c > $OUT_FOLDER/$SAMPLE.mapped_UMIs_per_barcode.txt

samtools view --threads $THREADS $BAM_FILE | awk '/xf:i:25/ && /CB:Z:/ {for(i=1;i<=NF;i++) if($i ~ /^CB:Z:/) {split($i,a,":"); print a[3]}}' | sort --parallel=$THREADS | uniq -c > $OUT_FOLDER/$SAMPLE.mapped_UMIs_per_barcode.txt

# get mapped reads counts
echo "Getting mapped reads counts"
#samtools view --threads $THREADS $BAM_FILE  | awk '($5 == 255 || $1 ~ /^@/)' | sed 's/.*CB:Z:\([ACGT]*\).*/\1/' | sort | uniq -c > $OUT_FOLDER/$SAMPLE.mapped_subset.txt
samtools view --threads $THREADS $BAM_FILE | awk '($5 == 255 || $1 ~ /^@/) && /CB:Z:/ {for(i=1;i<=NF;i++) if($i ~ /^CB:Z:/) {split($i,a,":"); print a[3]}}' | sort --parallel=$THREADS | uniq -c > $OUT_FOLDER/$SAMPLE.mapped_subset.txt

echo "Done!"