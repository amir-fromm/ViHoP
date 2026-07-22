#!/usr/bin/env bash
# Combine BLAST (+ optional LAST) host annotations into summarized tables.
# When LAST summary is missing / run_lastal is off, produce a BLAST-only summary
# with the column names expected by bin/blast_viruses_summary.py.
#
# Usage: run_summarize_annotations.sh BLAST_SUMMARY.csv CELLS.tsv OUT_FULL.csv \
#          OUT_SUM.csv LEVEL PERCENT [LAST_SUMMARY.csv] [LAST_COUNT]
set -euo pipefail

blast_summary="${1:?}"
cells_tsv="${2:?}"
out_full="${3:?}"
out_sum="${4:?}"
level="${5:-Class}"
percent="${6:-80}"
last_summary="${7:-}"
last_count="${8:-10}"

root="$(cd "$(dirname "$0")/../.." && pwd)"
mkdir -p "$(dirname "$out_full")" "$(dirname "$out_sum")"

# Effective last-read threshold (legacy doubles LAST_COUNT).
last_count_final=$(( last_count * 2 ))

if [[ -n "$last_summary" && -s "$last_summary" ]]; then
    (
        cd "$root"
        python bin/summarize_annotations.py \
            --last_path "$last_summary" \
            --blast_path "$blast_summary" \
            --cell_path "$cells_tsv" \
            --count_reads "$last_count_final" \
            --percent "$percent" \
            --level "$level"
        mv -f "full_annotation.${level}.csv" "$out_full"
        mv -f "summarized_annotation.${level}.csv" "$out_sum"
    )
else
    echo "[summarize] LAST summary absent; writing BLAST-only annotation"
    python - "$blast_summary" "$cells_tsv" "$out_full" "$out_sum" "$level" "$percent" <<'PY'
import sys
from pathlib import Path

import pandas as pd

blast_path, cells_path, out_full, out_sum, level, percent = sys.argv[1:7]
percent = float(percent) / 100.0

blast = pd.read_csv(blast_path)
cells = pd.read_table(cells_path, names=["Barcode", "Sample"])
cells["cell"] = cells["Barcode"] + "." + cells["Sample"]

if level not in blast.columns:
    raise SystemExit(f"ERROR: level column '{level}' missing from {blast_path}; columns={list(blast.columns)}")

blast = blast.rename(
    columns={
        level: "Blast Annotation",
        "length": "Blast max alignment length",
        "percent": "Blast percent alignment",
    }
)
blast[["Barcode", "Sample"]] = blast["cell"].str.split(".", n=1, expand=True)
blast["last Annotation"] = pd.NA
blast["last number of mapped reads"] = pd.NA
blast["last percent mapped reads"] = pd.NA

cols = [
    "cell",
    "Barcode",
    "Sample",
    "Blast Annotation",
    "Blast percent alignment",
    "Blast max alignment length",
    "last Annotation",
    "last number of mapped reads",
    "last percent mapped reads",
]
full = blast[cols].copy()
full = pd.merge(full, cells[["cell", "Barcode", "Sample"]], on=["cell", "Barcode", "Sample"], how="outer")
full = full.sort_values(
    ["cell", "Blast percent alignment"],
    ascending=[True, False],
)

# Keep top annotation per cell that passes chimeric threshold.
sum_df = full.dropna(subset=["Blast Annotation"]).copy()
sum_df = sum_df[sum_df["Blast percent alignment"] >= percent]
sum_df = sum_df.sort_values("Blast percent alignment", ascending=False).drop_duplicates("cell")

Path(out_full).parent.mkdir(parents=True, exist_ok=True)
full.to_csv(out_full, index=False)
sum_df.to_csv(out_sum, index=False)
print(f"[summarize] BLAST-only: {sum_df.shape[0]} cells annotated -> {out_sum}")
PY
fi

echo "[summarize] wrote $out_sum"
