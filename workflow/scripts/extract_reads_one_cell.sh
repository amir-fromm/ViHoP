#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "Usage: extract_reads_one_cell.sh --barcode BARCODE --sample SAMPLE --fastqs FASTQ_DIRS --outdir OUTDIR"
}

barcode=""
sample=""
fastqs=""
outdir=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --barcode) barcode="$2"; shift 2 ;;
        --sample) sample="$2"; shift 2 ;;
        --fastqs) fastqs="$2"; shift 2 ;;
        --outdir) outdir="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1"; usage; exit 1 ;;
    esac
done

if [[ -z "$barcode" || -z "$sample" || -z "$fastqs" || -z "$outdir" ]]; then
    usage
    exit 1
fi

mkdir -p "$outdir"

barcode_clean="$(echo "$barcode" | tr -d '[:space:]')"
tmpdir="$(mktemp -d)"
raw_out="$outdir/reads.raw.fq.gz"

trap 'rm -rf "$tmpdir"' EXIT

: > "$tmpdir/all_reads.fq"

IFS=',' read -ra fastq_dirs <<< "$fastqs"

for fastq_dir in "${fastq_dirs[@]}"; do
    for r1 in "$fastq_dir"/*_R1_*.fastq.gz; do
        [[ -e "$r1" ]] || continue

        r2="${r1/_R1_/_R2_}"
        [[ -e "$r2" ]] || {
            echo "Warning: missing R2 for $r1" >&2
            continue
        }

        base="$(basename "$r1" .fastq.gz)"
        headers="$tmpdir/${base}.headers.txt"

        zgrep -B1 "$barcode_clean" "$r1" \
            | grep '^@' \
            | cut -d' ' -f1 \
            | tr -d '@' \
            > "$headers"

        if [[ ! -s "$headers" ]]; then
            echo "No reads found for $barcode_clean in $r1" >&2
            continue
        fi

        seqtk subseq "$r2" "$headers" >> "$tmpdir/all_reads.fq"
    done
done

if [[ ! -s "$tmpdir/all_reads.fq" ]]; then
    echo "No reads extracted for $barcode_clean" >&2
    exit 1
fi

gzip -c "$tmpdir/all_reads.fq" > "$raw_out"

cat > "$outdir/metadata.yaml" <<EOF
cell: ${barcode_clean}.${sample}
barcode: ${barcode_clean}
sample: ${sample}
files:
  raw_reads: reads.raw.fq.gz
EOF
