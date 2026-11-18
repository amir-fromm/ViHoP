#!/usr/bin/env bash
# Sanitize FASTA headers for Cell Ranger (legacy RUN_DB_EDIT clean steps).
# Usage: run_clean_fasta.sh INPUT.fasta OUTPUT.fasta
set -euo pipefail

infile="${1:?}"
outfile="${2:?}"

mkdir -p "$(dirname "$outfile")"
tmp="$(mktemp "${TMPDIR:-/tmp}/clean_fasta.XXXXXX")"

# Replace characters Cell Ranger / GTF tooling dislike, then linearize sequences.
tr '|=()' '#' < "$infile" > "$tmp.clean"
awk '/^>/ {printf("\n%s\n",$0); next} {printf("%s",$0)} END {printf("\n")}' \
    < "$tmp.clean" > "$tmp.final"
# Drop leading blank line from awk.
sed '1{/^[[:space:]]*$/d;}' "$tmp.final" > "$outfile"

rm -f "$tmp" "$tmp.clean" "$tmp.final"
echo "[clean] wrote $outfile ($(grep -c '^>' "$outfile" || true) seqs)"
