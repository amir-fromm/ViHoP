#!/usr/bin/env bash
# Virus marker DIAMOND blastx, then general DIAMOND on those candidates.
# Same logic as bin/lastal_search_viruses.sh: drop queries whose best general
# hit is non-viral, then pick the best bitscore across marker + general hits.
# Snakemake owns the job. Do not pass --cluster lsf to diamond_to_general.sh.
#
# Usage: run_virus_diamond.sh TRANSCRIPTS.fasta MARKER.dmnd GENERAL.dmnd OUT_SUMMARY.tsv \
#          THREADS WORKDIR CATALOG NODES MERGED [DIAMOND_BIN]
set -euo pipefail

transcripts="${1:?}"
marker_db="${2:?}"
diamond_db="${3:?}"
out_summary="${4:?}"
threads="${5:-16}"
workdir="${6:-results/virus}"
catalog="${7:?}"
nodes="${8:?}"
merged="${9:?}"
diamond_bin="${10:-diamond}"

root="$(cd "$(dirname "$0")/../.." && pwd)"
mkdir -p "$workdir" "$(dirname "$out_summary")"

if [[ ! -s "$transcripts" ]] || ! grep -q '^>' "$transcripts"; then
    : > "$out_summary"
    echo "[virus] empty batch $transcripts; empty $out_summary"
    exit 0
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

marker_tsv="$workdir/marker_diamond.tsv"
marker_tmp="$workdir/marker_tmp"
mkdir -p "$marker_tmp"
echo "[virus] diamond blastx vs marker DB: $marker_db"
"$diamond_bin" blastx \
    --db "$marker_db" \
    --query "$transcripts" \
    --out "$marker_tsv" \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
    --threads "$threads" \
    --max-target-seqs 1 \
    --max-hsps 1 \
    --evalue 0.001 \
    --tmpdir "$marker_tmp" \
    --very-sensitive

awk '!/^#/ {print $1}' "$marker_tsv" | sort -u > "$workdir/marker_candidates.txt"
seqtk subseq "$transcripts" "$workdir/marker_candidates.txt" > "$workdir/marker_candidates.fasta"

n_cand=$(grep -c '^>' "$workdir/marker_candidates.fasta" || true)
echo "[virus] marker candidates: $n_cand"
if [[ "${n_cand:-0}" -eq 0 ]]; then
    : > "$out_summary"
    echo "[virus] no marker hits; empty $out_summary"
    exit 0
fi

diamond_dir="$workdir/diamond_general"
# A leftover alignment would make diamond_to_general.sh skip the search.
rm -rf "$diamond_dir"

echo "[virus] diamond blastx vs general DB: $diamond_db"
bash "$root/workflow/scripts/diamond_to_general.sh" \
    --cluster none \
    --input-fasta "$workdir/marker_candidates.fasta" \
    --output-dir "$diamond_dir" \
    --database "$diamond_db" \
    --catalog "$catalog" \
    --nodes "$nodes" \
    --merged "$merged" \
    --diamond "$diamond_bin" \
    --threads "$threads" \
    --sensitivity very-sensitive

general_aln="$diamond_dir/diamond_candidates_all.tsv"
false_positives="$diamond_dir/false_positives.list"
: > "$false_positives"
if [[ -s "$diamond_dir/edited_diamond_all.tsv" ]]; then
    # Queries whose best general hit is not Viruses (same as lastal_search_viruses.sh).
    awk -F '\t' '
        !/^#/ && ($14 == "Bacteria" || $14 == "Eukaryotes" || $14 == "Others" || $14 == "Phages") {
            if (!seen[$1]++) print $1
        }
    ' "$diamond_dir/edited_diamond_all.tsv" > "$false_positives"
fi

# Merge general + marker; best bitscore across both wins (marker subjects not in
# the RefSeq catalog are labeled Viruses by summarize_diamonds.py).
combined="$workdir/diamond_combined.tsv"
{
    if [[ -s "$general_aln" ]]; then
        grep -v '^#' "$general_aln" || true
    fi
    if [[ -s "$marker_tsv" ]]; then
        grep -v '^#' "$marker_tsv" || true
    fi
} > "$combined"

if [[ ! -s "$combined" ]]; then
    : > "$out_summary"
    echo "[virus] no combined DIAMOND hits; empty $out_summary"
    exit 0
fi

edited_combined="$workdir/edited_diamond_combined.tsv"
line_count=$(wc -l < "$combined")
echo "[virus] best-hit summary on marker+general ($line_count rows)"
python "$root/workflow/scripts/summarize_diamonds.py" \
    --alignment "$combined" \
    --false_positives "$false_positives" \
    --output_path "$edited_combined" \
    --catalog "$catalog" \
    --nodes "$nodes" \
    --merged "$merged" \
    --size_table "$line_count"

# 12 BLAST columns, virus-division best hits only. Pairing reads this file.
awk -F '\t' '!/^#/ && $14 == "Viruses" {
    OFS = "\t"
    print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12
}' "$edited_combined" > "$out_summary"

echo "[virus] wrote $out_summary ($(wc -l < "$out_summary") hits)"
