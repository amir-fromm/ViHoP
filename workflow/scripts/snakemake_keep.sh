#!/usr/bin/env bash
# Sidecar copies for expensive outputs. Snakemake deletes rule outputs before
# a re-run, so an in-place "if -s output" check always fails. Files under a
# sibling _keep/ directory are not rule outputs and survive that wipe.
#
# Usage:
#   snakemake_keep.sh restore KEEP_DIR FILE [FILE ...]
#     Copy each FILE from KEEP_DIR/$(basename FILE) when FILE is missing/empty
#     and the keep copy is non-empty.
#   snakemake_keep.sh save KEEP_DIR FILE [FILE ...]
#     Copy each non-empty FILE into KEEP_DIR.
#   snakemake_keep.sh all_present FILE [FILE ...]
#     Exit 0 if every FILE is non-empty, else 1.
set -euo pipefail

cmd="${1:?}"
shift

case "$cmd" in
  restore)
    keep_dir="${1:?}"
    shift
    mkdir -p "$keep_dir"
    for f in "$@"; do
      base="$(basename "$f")"
      keep="$keep_dir/$base"
      if [[ ! -s "$f" && -s "$keep" ]]; then
        mkdir -p "$(dirname "$f")"
        cp -a "$keep" "$f"
        echo "[keep] restored $f from $keep"
      fi
    done
    ;;
  save)
    keep_dir="${1:?}"
    shift
    mkdir -p "$keep_dir"
    for f in "$@"; do
      if [[ -s "$f" ]]; then
        cp -a "$f" "$keep_dir/$(basename "$f")"
      fi
    done
    ;;
  all_present)
    for f in "$@"; do
      if [[ ! -s "$f" ]]; then
        exit 1
      fi
    done
    exit 0
    ;;
  *)
    echo "Usage: $0 restore|save|all_present ..." >&2
    exit 2
    ;;
esac
