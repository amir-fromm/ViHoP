#!/usr/bin/env bash
# Nature software-checklist demo: host–virus pairing on two real Gossenkollesee cells.
# Inputs under demo/input/ include combined FASTQs from the manuscript run plus
# the lastal hits / host annotations needed for the pairing stage.
# Does not require Cell Ranger, SPAdes, or a cluster.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DEMO="$ROOT/demo"
OUT="${TMPDIR:-/tmp}/vihop_demo_out_$$"
mkdir -p "$OUT"
cleanup() { rm -rf "$OUT"; }
trap cleanup EXIT

echo "[demo] cells:"
cat "$DEMO/input/demo_cells.txt"
echo "[demo] FASTQs:"
find "$DEMO/input/cells" -name '*.combined.fq.gz' -exec ls -lh {} \;

python3 "$ROOT/bin/blast_viruses_summary.py" \
  --method blast \
  --output_path "$OUT" \
  --taxa_file "$DEMO/input/virus_taxonomy.csv" \
  --blast_file "$DEMO/input/virus_lastal.tsv" \
  --annotation_file "$DEMO/input/blast_summary.Order.csv" \
  --identity 70 \
  --length 30 \
  --output_name summary_blast_viruses.Order.tsv \
  --taxa Order \
  --threshold 100

echo "[demo] comparing to expected/"
diff -u "$DEMO/expected/summary_blast_viruses.Order.tsv" \
        "$OUT/summary_blast_viruses.Order.tsv"
diff -u "$DEMO/expected/virus_contigs.txt" \
        "$OUT/virus_contigs.txt"

echo "[demo] PASS: pairing summary matches expected/ for the two Gossenkollesee demo cells"
