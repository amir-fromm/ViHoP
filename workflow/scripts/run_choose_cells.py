#!/usr/bin/env python3
"""Build candidate_cells.tsv from Cell Ranger filtered matrices.

Default: every barcode in each sample's filtered_feature_bc_matrix.
Optional overrides: --min-umi (UMIs/cell), --min-genes (genes detected/cell).

Output (no header): barcode<TAB>sample
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path


def _open_text(path: Path):
    if path.suffix == ".gz" or str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open()


def read_barcodes(matrix_dir: Path) -> list[str]:
    for name in ("barcodes.tsv.gz", "barcodes.tsv"):
        path = matrix_dir / name
        if path.exists():
            with _open_text(path) as handle:
                return [line.strip().split("\t")[0] for line in handle if line.strip()]
    raise FileNotFoundError(f"No barcodes.tsv(.gz) in {matrix_dir}")


def strip_gem_suffix(barcode: str) -> str:
    # Cell Ranger barcodes look like ACGTACGTACGTACGT-1
    if barcode.endswith("-1"):
        return barcode[: -2]
    return barcode.split("-")[0] if "-" in barcode else barcode


def load_counts(matrix_dir: Path):
    """Return (barcodes, sparse matrix genes x cells) from a 10x matrix dir."""
    try:
        from scipy.io import mmread
    except ImportError as exc:
        raise SystemExit(
            "ERROR: scipy is required when min_umi / min_genes filters are set. "
            "Install scipy in the vihop env, or leave filters null to keep all "
            "filtered Cell Ranger barcodes."
        ) from exc

    barcodes = read_barcodes(matrix_dir)
    mtx_path = matrix_dir / "matrix.mtx.gz"
    if not mtx_path.exists():
        mtx_path = matrix_dir / "matrix.mtx"
    if not mtx_path.exists():
        raise FileNotFoundError(f"No matrix.mtx(.gz) in {matrix_dir}")
    mat = mmread(str(mtx_path)).tocsc()
    if mat.shape[1] != len(barcodes):
        raise ValueError(
            f"Matrix columns ({mat.shape[1]}) != barcodes ({len(barcodes)}) in {matrix_dir}"
        )
    return barcodes, mat


def select_cells(
    matrix_dir: Path,
    min_umi: int | None,
    min_genes: int | None,
) -> list[str]:
    """Return GEM-stripped barcodes passing filters (or all if filters are None)."""
    need_matrix = (min_umi is not None and min_umi > 0) or (
        min_genes is not None and min_genes > 0
    )
    if not need_matrix:
        return [strip_gem_suffix(bc) for bc in read_barcodes(matrix_dir)]

    import numpy as np

    barcodes, mat = load_counts(matrix_dir)
    # mat: genes x cells
    keep = np.ones(mat.shape[1], dtype=bool)
    if min_umi is not None and min_umi > 0:
        umi = np.asarray(mat.sum(axis=0)).ravel()
        keep &= umi >= min_umi
    if min_genes is not None and min_genes > 0:
        n_genes = np.asarray((mat > 0).sum(axis=0)).ravel()
        keep &= n_genes >= min_genes
    return [strip_gem_suffix(bc) for bc, ok in zip(barcodes, keep) if ok]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--bam-root",
        required=True,
        help="Parent of <sample>/outs/filtered_feature_bc_matrix",
    )
    p.add_argument(
        "--samples",
        required=True,
        help="Comma-separated sample IDs",
    )
    p.add_argument(
        "--output",
        required=True,
        help="candidate_cells.tsv path (barcode<TAB>sample)",
    )
    p.add_argument(
        "--min-umi",
        default=None,
        type=int,
        help="Min UMIs per cell (omit / 0 = keep all filtered cells)",
    )
    p.add_argument(
        "--min-genes",
        default=None,
        type=int,
        help="Min genes detected per cell (omit / 0 = keep all filtered cells)",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    samples = [s.strip() for s in args.samples.split(",") if s.strip()]
    if not samples:
        raise SystemExit("ERROR: empty --samples")

    min_umi = args.min_umi if args.min_umi and args.min_umi > 0 else None
    min_genes = args.min_genes if args.min_genes and args.min_genes > 0 else None

    bam_root = Path(args.bam_root)
    rows: list[tuple[str, str]] = []
    for sample in samples:
        matrix_dir = bam_root / sample / "outs" / "filtered_feature_bc_matrix"
        if not matrix_dir.is_dir():
            raise FileNotFoundError(f"Missing filtered matrix: {matrix_dir}")
        cells = select_cells(matrix_dir, min_umi, min_genes)
        print(
            f"[choose_cells] {sample}: {len(cells)} cells "
            f"(min_umi={min_umi or 'all'}, min_genes={min_genes or 'all'})"
        )
        for bc in cells:
            rows.append((bc, sample))

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as handle:
        for bc, sample in rows:
            handle.write(f"{bc}\t{sample}\n")
    print(f"[choose_cells] wrote {len(rows)} cells → {out}")


if __name__ == "__main__":
    main()
