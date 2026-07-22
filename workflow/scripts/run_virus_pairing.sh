#!/usr/bin/env bash
# Wrap bin/blast_viruses_summary.py for host–virus pairing.
# Usage: run_virus_pairing.sh VIRUS_HITS.tsv ANNOTATION.csv TAXA.csv OUT_DIR OUT_NAME \
#          IDENTITY LENGTH
set -euo pipefail

virus_hits="${1:?}"
annotation="${2:?}"
taxa="${3:?}"
out_dir="${4:?}"
out_name="${5:?}"
identity="${6:-70}"
length="${7:-30}"

root="$(cd "$(dirname "$0")/../.." && pwd)"
mkdir -p "$out_dir"

if [[ ! -s "$virus_hits" ]]; then
    echo "[virus_pair] empty virus hits; writing empty $out_dir/$out_name"
    : > "$out_dir/$out_name"
    : > "$out_dir/virus_contigs.txt"
    exit 0
fi

python "$root/bin/blast_viruses_summary.py" \
    --output_path "$out_dir" \
    --output_name "$out_name" \
    --taxa_file "$taxa" \
    --blast_file "$virus_hits" \
    --annotation_file "$annotation" \
    --identity "$identity" \
    --length "$length"

echo "[virus_pair] wrote $out_dir/$out_name"
