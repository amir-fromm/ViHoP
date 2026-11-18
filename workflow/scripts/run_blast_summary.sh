#!/usr/bin/env bash
# Filter BLAST TSVs and summarize with bin/blast_sum.py (TSV path; no parquet deps).
# Usage: run_blast_summary.sh PR2.tsv METAPR2.tsv OUT_SUMMARY.csv REPLACEMENT.tsv \
#          IDENTITY LENGTH LEVEL [PREFILTER_PIDENT] [PREFILTER_LENGTH] [FORCE_ANNOTATION]
# FORCE_ANNOTATION defaults to true: a real genus beats a placeholder even at a
# lower bitscore, and placeholder contigs are dropped when the cell has a real genus.
set -euo pipefail

pr2_tsv="${1:?}"
metapr2_tsv="${2:?}"
out_summary="${3:?}"
replacement="${4:?}"
identity="${5:-99}"
length="${6:-300}"
level="${7:-Class}"
pre_pident="${8:-90}"
pre_length="${9:-100}"
force_annotation="${10:-true}"

root="$(cd "$(dirname "$0")/../.." && pwd)"
workdir="$(pwd)"
mkdir -p "$(dirname "$out_summary")"

# Filtered tables stay in this work directory. A named annotation run
# (blast_summary.Genus.100.csv) writes its own filtered tables so the
# unsuffixed 300 bp tables are left in place.
summary_base="$(basename "$out_summary")"
if [[ "$summary_base" == "blast_summary.${level}.csv" ]]; then
    filt_pr2="$workdir/results/host/filtered.rRNA.PR2.tsv"
    filt_meta="$workdir/results/host/filtered.rRNA.metaPR2.tsv"
else
    tag="${summary_base#blast_summary.}"
    tag="${tag%.csv}"
    filt_pr2="$workdir/results/host/filtered.rRNA.PR2.${tag}.tsv"
    filt_meta="$workdir/results/host/filtered.rRNA.metaPR2.${tag}.tsv"
fi
if [[ "$replacement" != /* ]]; then
    if [[ -f "$workdir/$replacement" ]]; then
        replacement="$workdir/$replacement"
    else
        replacement="$root/$replacement"
    fi
fi

# Apply summary cutoffs while streaming, strip barrnap :start-end(strand)
# from qseqid, and keep nuclear/nucleomorph 18S only for PR2 (drop plastid
# 16S and other non-18S PR2 subjects). metaPR2 has no organelle field here.
echo "Filtering PR2 (identity >= ${identity}, length >= ${length}, 18S_rRNA, not plastid, force_annotation=${force_annotation})"
awk -v p="$identity" -v l="$length" '
  BEGIN { FS = OFS = "\t" }
  $3 >= p && $4 >= l {
    n = split($2, a, "|")
    if (n >= 3 && a[2] == "18S_rRNA" && a[3] != "plastid") {
      sub(/:[0-9]+-[0-9]+\([+-]\)$/, "", $1)
      print
    }
  }
' "$pr2_tsv" > "$filt_pr2"

echo "Filtering metaPR2 (identity >= ${identity}, length >= ${length})"
awk -v p="$identity" -v l="$length" '
  BEGIN { FS = OFS = "\t" }
  $3 >= p && $4 >= l {
    sub(/:[0-9]+-[0-9]+\([+-]\)$/, "", $1)
    print
  }
' "$metapr2_tsv" > "$filt_meta"

mkdir -p "$(dirname "$out_summary")"
# blast_sum.py always writes blast_summary.<level>.csv into --path.
# Use a private directory so that name is not the unsuffixed summary.
tmp_sum="$(mktemp -d)"
trap 'rm -rf "$tmp_sum"' EXIT

python "$root/bin/blast_sum.py" \
    --path "$tmp_sum" \
    --metapr2_file "$filt_meta" \
    --pr2_file "$filt_pr2" \
    --replacement_file "$replacement" \
    --identity "$identity" \
    --length "$length" \
    --level "$level" \
    --force-annotation "$force_annotation"

mv -f "$tmp_sum/blast_summary.${level}.csv" "$out_summary"
echo "[blast_summary] wrote $out_summary ($(wc -l < "$out_summary") lines)"
