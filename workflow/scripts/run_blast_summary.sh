#!/usr/bin/env bash
# Filter BLAST TSVs and summarize with bin/blast_sum.py (TSV path; no parquet deps).
# Usage: run_blast_summary.sh PR2.tsv METAPR2.tsv OUT_SUMMARY.csv REPLACEMENT.tsv \
#          IDENTITY LENGTH LEVEL [PREFILTER_PIDENT] [PREFILTER_LENGTH]
set -euo pipefail

pr2_tsv="${1:?}"
metapr2_tsv="${2:?}"
out_summary="${3:?}"
replacement="${4:?}"
identity="${5:-98}"
length="${6:-300}"
level="${7:-Class}"
pre_pident="${8:-90}"
pre_length="${9:-100}"

root="$(cd "$(dirname "$0")/../.." && pwd)"
mkdir -p "$(dirname "$out_summary")" results/host

filt_pr2="results/host/filtered.rRNA.PR2.tsv"
filt_meta="results/host/filtered.rRNA.metaPR2.tsv"

# Prefilter, and strip barrnap :start-end(strand) from qseqid so blast_sum.py
# can parse cell = barcode.sample from ..._i0_<cell>.
awk -v p="$pre_pident" -v l="$pre_length" '
  $3 >= p && $4 >= l {
    sub(/:[0-9]+-[0-9]+\([+-]\)$/, "", $1)
    print
  }
' OFS='\t' "$pr2_tsv" > "$filt_pr2"

awk -v p="$pre_pident" -v l="$pre_length" '
  $3 >= p && $4 >= l {
    sub(/:[0-9]+-[0-9]+\([+-]\)$/, "", $1)
    print
  }
' OFS='\t' "$metapr2_tsv" > "$filt_meta"

tmp_summary="$root/blast_summary.${level}.csv"
rm -f "$tmp_summary"

# blast_sum.py joins --path + file (file must start with /) and writes
# blast_summary.{level}.csv under --path.
python "$root/bin/blast_sum.py" \
    --path "$root" \
    --metapr2_file "/${filt_meta}" \
    --pr2_file "/${filt_pr2}" \
    --replacement_file "$replacement" \
    --identity "$identity" \
    --length "$length" \
    --level "$level"

mv -f "$tmp_summary" "$out_summary"
echo "[blast_summary] wrote $out_summary ($(wc -l < "$out_summary") lines)"
