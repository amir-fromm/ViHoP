#!/usr/bin/env bash
# Run cellranger count (no nested bsub). Supports one or more FASTQ dirs.
# Usage:
#   run_cellranger_count.sh SAMPLE TRANSCRIPTOME FASTQS_CSV OUTDIR THREADS MEM_MB \
#       [CREATE_BAM] [CHEMISTRY] [EXPECT_CELLS] [MODULE] [LOCALMEM_GB] [LOCALVMEM_GB] \
#       [OVERRIDES_JSON] [FORCE_CLEAN]
#
# FASTQS_CSV: single path or comma-separated paths (Cell Ranger native).
# EXPECT_CELLS: empty / "false" / "null" → omit --expect-cells.
# MODULE: e.g. cellranger/10.1.0 (empty → try CellRanger/9.0.0 then CellRanger).
# LOCALMEM_GB / LOCALVMEM_GB: empty → derive localmem from MEM_MB; omit localvmem.
# OVERRIDES_JSON: path to Martian --overrides JSON, or empty to omit.
# FORCE_CLEAN: true → rm -rf sample pipestance before run; otherwise resume
# only when _invocation is present. Snakemake creates the sample directory
# before the job because the BAM and barcodes live inside it. That directory
# is not a pipestance, so it is removed and count starts fresh.
set -euo pipefail

sample="${1:?}"
transcriptome="${2:?}"
fastqs="${3:?}"
outdir="${4:?}"
threads="${5:-16}"
mem_mb="${6:-64000}"
create_bam="${7:-true}"
chemistry="${8:-auto}"
expect_cells="${9:-}"
module_name="${10:-}"
localmem_gb="${11:-}"
localvmem_gb="${12:-}"
overrides_json="${13:-}"
force_clean="${14:-false}"

ulimit -v unlimited || true

load_cellranger() {
    if command -v cellranger >/dev/null 2>&1; then
        return 0
    fi
    if ! command -v module >/dev/null 2>&1; then
        return 1
    fi
    module purge || true
    if [[ -n "$module_name" && "$module_name" != "null" && "$module_name" != "false" ]]; then
        module load "$module_name" || true
    else
        module load CellRanger/9.0.0 || module load CellRanger || module load cellranger/10.1.0 || true
    fi
    command -v cellranger >/dev/null 2>&1
}

if ! load_cellranger; then
    echo "ERROR: cellranger not found (module=${module_name:-default})" >&2
    exit 1
fi

# Absolute transcriptome path — we cd into outdir before invoking cellranger.
if [[ ! -d "$transcriptome" ]]; then
    echo "ERROR: transcriptome missing: $transcriptome" >&2
    exit 1
fi
transcriptome="$(cd "$transcriptome" && pwd)"

IFS=',' read -r -a fq_arr <<< "$fastqs"
abs_fastqs=()
for d in "${fq_arr[@]}"; do
    if [[ ! -d "$d" ]]; then
        echo "ERROR: FASTQ directory missing: $d" >&2
        exit 1
    fi
    abs="$(cd "$d" && pwd)"
    n_fq=$(find "$abs" -maxdepth 1 -type f \( -name '*.fastq.gz' -o -name '*.fq.gz' \) | wc -l | tr -d ' ')
    if [[ "$n_fq" -eq 0 ]]; then
        echo "ERROR: no *.fastq.gz / *.fq.gz in $abs" >&2
        echo "Directory listing:" >&2
        ls -la "$abs" >&2 || true
        exit 1
    fi
    echo "[cellranger] $abs ($n_fq fastq files)"
    abs_fastqs+=("$abs")
done
fastqs="$(IFS=','; echo "${abs_fastqs[*]}")"

mkdir -p "$outdir"
outdir="$(cd "$outdir" && pwd)"
sample_dir="$outdir/$sample"
bam="$sample_dir/outs/possorted_genome_bam.bam"
matrix_dir="$sample_dir/outs/filtered_feature_bc_matrix"
barcodes="$matrix_dir/barcodes.tsv.gz"

if [[ -f "$barcodes" && ( -f "$bam" || "$create_bam" == "false" ) ]]; then
    echo "[cellranger] existing outs for $sample — skipping"
    exit 0
fi

if [[ "$force_clean" == "true" ]]; then
    if [[ -d "$sample_dir" ]]; then
        echo "[cellranger] force_clean: removing $sample_dir"
        rm -rf "$sample_dir"
    fi
elif [[ -d "$sample_dir" && ! -f "$sample_dir/_invocation" ]]; then
    echo "[cellranger] $sample_dir is not a pipestance — removing it so count can start"
    rm -rf "$sample_dir"
elif [[ -d "$sample_dir" ]]; then
    # Resume: drop stale Martian lock only if no live cellranger process.
    if [[ -e "$sample_dir/_lock" ]]; then
        if pgrep -f "cellranger.*--id[= ]${sample}" >/dev/null 2>&1; then
            echo "ERROR: $sample_dir/_lock present and cellranger still running" >&2
            exit 1
        fi
        echo "[cellranger] removing stale _lock for resume"
        rm -f "$sample_dir/_lock"
    fi
    echo "[cellranger] resuming existing pipestance $sample_dir"
fi

export TENX_DISABLE_UI="${TENX_DISABLE_UI:-1}"
export MRO_DISK_SPACE_CHECK="${MRO_DISK_SPACE_CHECK:-disable}"

if [[ -z "$localmem_gb" || "$localmem_gb" == "null" || "$localmem_gb" == "false" ]]; then
    localmem_gb=$(( mem_mb / 1000 ))
fi
if [[ "$localmem_gb" -lt 1 ]]; then
    localmem_gb=1
fi

cmd=(
    cellranger count
    --id="$sample"
    --transcriptome="$transcriptome"
    --fastqs="$fastqs"
    --create-bam="$create_bam"
    --nosecondary
    --localcores="$threads"
    --localmem="$localmem_gb"
)
if [[ -n "$localvmem_gb" && "$localvmem_gb" != "null" && "$localvmem_gb" != "false" ]]; then
    cmd+=(--localvmem="$localvmem_gb")
fi
if [[ -n "$overrides_json" && "$overrides_json" != "null" && "$overrides_json" != "false" ]]; then
    if [[ ! -f "$overrides_json" ]]; then
        echo "ERROR: overrides JSON missing: $overrides_json" >&2
        exit 1
    fi
    overrides_abs="$(cd "$(dirname "$overrides_json")" && pwd)/$(basename "$overrides_json")"
    cmd+=(--overrides="$overrides_abs")
fi
if [[ -n "$chemistry" && "$chemistry" != "null" && "$chemistry" != "false" ]]; then
    cmd+=(--chemistry="$chemistry")
fi
if [[ -n "$expect_cells" && "$expect_cells" != "null" && "$expect_cells" != "false" ]]; then
    cmd+=(--expect-cells="$expect_cells")
fi

n_fq=${#fq_arr[@]}
echo "[cellranger] sample=$sample fastq_dirs=$n_fq transcriptome=$transcriptome"
echo "[cellranger] fastqs=$fastqs"
echo "[cellranger] module=${module_name:-default} localcores=$threads localmem=$localmem_gb localvmem=${localvmem_gb:-omit}"
echo "[cellranger] overrides=${overrides_json:-omit} force_clean=$force_clean"
echo "[cellranger] create_bam=$create_bam chemistry=${chemistry:-omit} expect_cells=${expect_cells:-omit}"
echo "[cellranger] ulimit -v=$(ulimit -v)"

(
    cd "$outdir"
    "${cmd[@]}"
)

bam="$sample_dir/outs/possorted_genome_bam.bam"
matrix_dir="$sample_dir/outs/filtered_feature_bc_matrix"
barcodes="$matrix_dir/barcodes.tsv.gz"
h5="$sample_dir/outs/filtered_feature_bc_matrix.h5"

publish_filtered_mex() {
    # Prefer a real Martian stage directory (ignore Mac SMB XSym text "symlinks").
    local stage_mex="" cand
    while IFS= read -r cand; do
        if [[ -d "$cand" && -f "$cand/barcodes.tsv.gz" ]]; then
            stage_mex="$cand"
            break
        fi
    done < <(find "$sample_dir/SC_RNA_COUNTER_CS" -name filtered_matrices_mex 2>/dev/null || true)

    if [[ -n "$stage_mex" ]]; then
        echo "[cellranger] publishing filtered MEX from Martian stage → outs/"
        mkdir -p "$matrix_dir"
        cp -a "$stage_mex/barcodes.tsv.gz" "$stage_mex/features.tsv.gz" "$stage_mex/matrix.mtx.gz" "$matrix_dir/"
        return 0
    fi

    # Last resort: barcodes from filtered H5 (features/matrix may already be in outs/).
    if [[ -f "$h5" && ! -f "$barcodes" ]]; then
        echo "[cellranger] restoring barcodes.tsv.gz from $h5"
        mkdir -p "$matrix_dir"
        python3 - "$h5" "$barcodes" <<'PY'
import gzip, sys
from pathlib import Path
h5, out = Path(sys.argv[1]), Path(sys.argv[2])
import h5py
with h5py.File(h5, "r") as f:
    bcs = f["matrix"]["barcodes"][:]
with gzip.open(out, "wt") as o:
    for b in bcs:
        o.write((b.decode() if isinstance(b, bytes) else str(b)) + "\n")
print(f"wrote {out} ({len(bcs)} barcodes)")
PY
    fi
}

# Cell Ranger sometimes leaves filtered MEX only under FILTER_BARCODES/.../filtered_matrices_mex
# (VDR publish to outs/ silently skipped). Rescue before declaring success.
if [[ ! -f "$barcodes" ]]; then
    publish_filtered_mex || true
fi

if [[ ! -f "$barcodes" ]]; then
    echo "ERROR: missing filtered barcodes at $barcodes" >&2
    exit 1
fi
if [[ "$create_bam" == "true" && ! -f "$bam" ]]; then
    echo "ERROR: missing BAM at $bam" >&2
    echo "ERROR: pipestance may be marked complete while outs/ BAM was deleted." >&2
    echo "ERROR: set params.cellranger.force_clean: true and re-run this sample." >&2
    exit 1
fi
echo "[cellranger] done $sample"
