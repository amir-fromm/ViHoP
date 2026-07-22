#!/usr/bin/env bash
# Virus LAST search on combined transcripts (no nested bsub; Snakemake owns the job).
# Mirrors bin/lastal_search_viruses.sh logic against results/ paths.
#
# Usage: run_virus_lastal.sh TRANSCRIPTS.fasta VIRUS_LAST_DB GENERAL_LAST_DB OUT_SUMMARY.tsv \
#          THREADS [WORKDIR]
set -euo pipefail

transcripts="${1:?}"
virus_db="${2:?}"
general_db="${3:?}"
out_summary="${4:?}"
threads="${5:-8}"
workdir="${6:-results/virus}"

root="$(cd "$(dirname "$0")/../.." && pwd)"
mkdir -p "$workdir"

if ! command -v lastal >/dev/null 2>&1; then
    if command -v module >/dev/null 2>&1; then
        module purge || true
        module load last || module load LAST || true
    fi
fi
if ! command -v lastal >/dev/null 2>&1; then
    echo "ERROR: lastal not found" >&2
    exit 1
fi
if ! command -v seqtk >/dev/null 2>&1; then
    if command -v module >/dev/null 2>&1; then
        module load seqtk || true
    fi
fi
if ! command -v seqtk >/dev/null 2>&1; then
    echo "ERROR: seqtk not found" >&2
    exit 1
fi

virus_last="$workdir/last_viruses.last"
candidates_txt="$workdir/last_viruses_candidates.txt"
candidates_fa="$workdir/last_viruses_candidates.fasta"
general_last="$workdir/last_candidates_all.last"
edited_all="$workdir/edited_last_all.tsv"
false_pos="$workdir/false_positives.list"
combined_last="$workdir/last_combined.last"
edited_combined="$workdir/edited_last_combined.tsv"

echo "[virus] lastal vs virus DB: $virus_db"
lastal -m 100 -f BlastTab -P "$threads" -i14G -u 2 -v -Q 0 -F 15 \
    "$virus_db" "$transcripts" > "$virus_last"

awk '!/^#/ {print $1}' "$virus_last" | sort -u > "$candidates_txt"
seqtk subseq "$transcripts" "$candidates_txt" > "$candidates_fa"

n_cand=$(grep -c '^>' "$candidates_fa" || true)
echo "[virus] candidate contigs: $n_cand"
if [[ "${n_cand:-0}" -eq 0 ]]; then
    : > "$out_summary"
    echo "[virus] no candidates; empty $out_summary"
    exit 0
fi

echo "[virus] lastal vs general DB: $general_db"
lastal -m 100 -f BlastTab -P "$threads" -i32G -u 2 -v -Q 0 -F 15 \
    "$general_db" "$candidates_fa" > "$general_last"

len=$(wc -l < "$general_last")
python "$root/bin/summarize_lastal_one_large.py" \
    --chunks 100000 \
    --size_table "$len" \
    --output_path "$edited_all" \
    --last_path "$general_last"

awk -F'\t' '
  !/^#/ && ($14 == "Bacteria" || $14 == "Eukaryotes" || $14 == "Others" || $14 == "Phages") {
    if (!seen[$1]++) print $1
  }
' "$edited_all" > "$false_pos"

{ grep -v '^#' "$general_last" || true
  grep -v '^#' "$virus_last" || true
} > "$combined_last"

len=$(wc -l < "$combined_last")
python "$root/bin/summarize_lastal_one_large.py" \
    --chunks 100000 \
    --size_table "$len" \
    --false_positives "$false_pos" \
    --output_path "$edited_combined" \
    --last_path "$combined_last"

awk -F'\t' '!/^#/ && $14 == "Viruses" {
    OFS="\t"
    print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12
}' "$edited_combined" > "$out_summary"

echo "[virus] wrote $out_summary ($(wc -l < "$out_summary") hits)"
