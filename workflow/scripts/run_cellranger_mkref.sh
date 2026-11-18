#!/usr/bin/env bash
# Run cellranger mkref without nested bsub (Snakemake owns the job).
# Usage: run_cellranger_mkref.sh REF_NAME FASTA GTF OUTDIR THREADS MEM_MB DONE_FLAG
set -euo pipefail

ref_name="${1:?}"
fasta="${2:?}"
gtf="${3:?}"
outdir="${4:?}"
threads="${5:-16}"
mem_mb="${6:-64000}"
done_flag="${7:?}"

# Guard against misconfigured outputs.cellranger_ref_done pointing at a real
# transcriptome reference.json (would overwrite the ref metadata with a date).
case "$(basename "$done_flag")" in
    reference.json|*.json)
        echo "ERROR: done_flag must be a .done marker, not $done_flag" >&2
        echo "       Set outputs.cellranger_ref_done to e.g. results/ref/NAME/mkref.done" >&2
        exit 1
        ;;
esac

if ! command -v cellranger >/dev/null 2>&1; then
    if command -v module >/dev/null 2>&1; then
        module purge || true
        module load CellRanger/9.0.0 || module load CellRanger || true
    fi
fi
if ! command -v cellranger >/dev/null 2>&1; then
    echo "ERROR: cellranger not found" >&2
    exit 1
fi

mem_gb=$(( mem_mb / 1000 ))
if [[ "$mem_gb" -lt 1 ]]; then
    mem_gb=1
fi

mkdir -p "$outdir"
ref_path="$outdir/$ref_name"
pipe_dir="$outdir/mkref_${ref_name}"
martian_ref="$pipe_dir/outs/reference"

ref_is_complete() {
    local d="$1"
    [[ -d "$d" && -f "$d/reference.json" ]]
}

# Reuse a finished reference (e.g. after Cell Ranger exited non-zero on cleanup).
if ref_is_complete "$ref_path"; then
    echo "[mkref] existing reference found at $ref_path — marking done"
    date > "$done_flag"
    exit 0
fi
if ref_is_complete "$martian_ref"; then
    echo "[mkref] promoting Martian outs/reference -> $ref_path"
    rm -rf "$ref_path"
    # Prefer a real directory for Cell Ranger count --transcriptome.
    mv "$martian_ref" "$ref_path"
    date > "$done_flag"
    exit 0
fi

if [[ ! -d "$pipe_dir" && -d "$ref_path" ]]; then
    echo "[mkref] removing incomplete $ref_path"
    rm -rf "$ref_path"
fi

export TENX_DISABLE_UI="${TENX_DISABLE_UI:-1}"
export MRO_DISK_SPACE_CHECK="${MRO_DISK_SPACE_CHECK:-disable}"

echo "[mkref] genome=$ref_name fasta=$fasta gtf=$gtf threads=$threads memgb=$mem_gb"
set +e
(
    cd "$outdir"
    cellranger mkref \
        --genome="$ref_name" \
        --fasta="$fasta" \
        --genes="$gtf" \
        --nthreads="$threads" \
        --memgb="$mem_gb"
)
rc=$?
set -e

# Cell Ranger 7–9 sometimes exits non-zero after a successful build
# (e.g. "ERROR: Directory not empty" while packing the .mri.tgz).
if ref_is_complete "$ref_path"; then
    echo "[mkref] reference OK at $ref_path (cellranger exit=$rc)"
    date > "$done_flag"
    exit 0
fi
if ref_is_complete "$martian_ref"; then
    echo "[mkref] promoting Martian outs/reference -> $ref_path (cellranger exit=$rc)"
    rm -rf "$ref_path"
    mv "$martian_ref" "$ref_path"
    date > "$done_flag"
    exit 0
fi

echo "ERROR: mkref failed (exit=$rc); no reference.json under:" >&2
echo "  $ref_path" >&2
echo "  $martian_ref" >&2
exit "$rc"
