#!/usr/bin/env bash
set -euo pipefail

FA=$1
OUT=$2
TMP=mm_tmp
THREADS=32

# Tunables (for "local duplication" detection):
MIN_ID=$3          # local identity threshold (e.g. 0.95–0.99)
COV=0.30             # fraction of query that must be covered by the alignment (local partial match)
COV_MODE=1           # 1 = query coverage (lets partial overlaps cluster)
MIN_ALN_LEN=200      # minimum aligned length in residues/bp (helps avoid tiny repeat matches)

mkdir -p "$OUT" "$TMP"

# Build DB
mmseqs createdb "$FA" "$OUT/db"

# Cluster (sensitive mode; local alignments are default)
mmseqs cluster \
  "$OUT/db" "$OUT/clu" "$TMP" \
  --min-seq-id "$MIN_ID" \
  -c "$COV" --cov-mode "$COV_MODE" \
  --min-aln-len "$MIN_ALN_LEN" \
  --threads "$THREADS"

# Export cluster membership (rep -> members)
mmseqs createtsv "$OUT/db" "$OUT/db" "$OUT/clu" "$OUT/clusters.tsv"

# Extract representative sequences (= deduplicated FASTA)
mmseqs result2repseq "$OUT/db" "$OUT/clu" "$OUT/repdb" --threads "$THREADS"
mmseqs result2flat  "$OUT/db" "$OUT/db" "$OUT/repdb" "$OUT/dedup_reps.fa" --use-fasta-header

echo "Done."
echo "Cluster table: $OUT/clusters.tsv"
echo "Deduplicated FASTA (representatives): $OUT/dedup_reps.fa"
