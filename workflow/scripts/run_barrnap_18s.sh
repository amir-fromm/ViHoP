#!/usr/bin/env bash
# Run barrnap on filtered transcripts and extract 18S / SSU sequences.
# Usage: run_barrnap_18s.sh INPUT.fasta RRNA.FA RRNA.GFF OUT_18S.FASTA THREADS [KINGDOM]
set -euo pipefail

infile="${1:?}"
rrna_fa="${2:?}"
rrna_gff="${3:?}"
out_18s="${4:?}"
threads="${5:-4}"
kingdom="${6:-fun}"

if ! command -v barrnap >/dev/null 2>&1; then
    if command -v module >/dev/null 2>&1; then
        module purge || true
        module load barrnap || true
    fi
fi

if ! command -v barrnap >/dev/null 2>&1; then
    echo "ERROR: barrnap not found in PATH or modules." >&2
    echo "Install into the vihop env, then re-run:" >&2
    echo "  conda install -c bioconda -c conda-forge barrnap" >&2
    exit 1
fi

if [[ ! -f "$infile" ]]; then
    echo "ERROR: missing input fasta: $infile" >&2
    exit 1
fi
if [[ ! -s "$infile" ]]; then
    echo "[barrnap] empty batch $infile; writing empty outputs"
    mkdir -p "$(dirname "$rrna_fa")" "$(dirname "$out_18s")"
    : > "$rrna_fa"
    : > "$rrna_gff"
    : > "$out_18s"
    exit 0
fi

mkdir -p "$(dirname "$rrna_fa")" "$(dirname "$out_18s")"

# barrnap writes barrnap.bed in the current directory and bedtools reads it
# back. Parallel LSF jobs share this work directory, so one job's bedtools
# sees another job's half-written BED ("less than 3 columns"). Run barrnap
# from a private directory. Paths must be absolute because barrnap chdirs.
abs_path() {
    local p="$1"
    if [[ "$p" == /* ]]; then
        printf '%s\n' "$p"
    else
        printf '%s\n' "$PWD/$p"
    fi
}
infile="$(abs_path "$infile")"
rrna_fa="$(abs_path "$rrna_fa")"
rrna_gff="$(abs_path "$rrna_gff")"
out_18s="$(abs_path "$out_18s")"

echo "[barrnap] version: $(barrnap --version 2>&1 || true)"
echo "[barrnap] in=$infile kingdom=$kingdom threads=$threads"

private="$(mktemp -d "${TMPDIR:-/tmp}/barrnap_job.XXXXXX")"
trap 'rm -rf "$private"' EXIT

# barrnap 1.x: kingdoms are bac|arc|fun (euk was renamed to fun).
# --rrna requests ribosomal RNA annotation; --outseq writes hit FASTA.
(
    cd "$private"
    barrnap \
        --kingdom "$kingdom" \
        --rrna \
        --threads "$threads" \
        --outseq "$rrna_fa" \
        "$infile" \
        > "$rrna_gff"
)

# barrnap 1.x puts rRNA type in GFF attributes (Name=18S_rRNA), not FASTA headers.
: > "$out_18s"
if [[ -s "$rrna_gff" && -s "$rrna_fa" ]]; then
    ids_file="$(mktemp "${TMPDIR:-/tmp}/barrnap_18s_ids.XXXXXX")"
    awk -F'\t' '
        $3 == "rRNA" && $9 ~ /Name=18S_rRNA|Alias=SSU_rRNA/ { ids[$1] = 1 }
        END { for (id in ids) print id }
    ' "$rrna_gff" | sort -u > "$ids_file"

    if [[ -s "$ids_file" ]]; then
        awk '
            NR == FNR { ids[$1] = 1; next }
            /^>/ {
                name = substr($0, 2)
                sub(/:.*/, "", name)
                keep = (name in ids)
            }
            keep { print }
        ' "$ids_file" "$rrna_fa" > "$out_18s"
    fi
    rm -f "$ids_file"
fi

echo "[barrnap] wrote $out_18s ($(grep -c '^>' "$out_18s" 2>/dev/null || echo 0) seqs)"
