#!/usr/bin/env python3
"""QC of host 18S annotation after the cell summary.

Inputs are the BLAST summary (all taxa per cell with alignment-length
percent), the pure-cell table, and the candidate cell list. Writes tables,
a text report, and SVG figures under an output directory.
"""

from __future__ import annotations

import argparse
import html
import math
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path

import pandas as pd


def load_candidates(path: Path) -> pd.DataFrame:
    cells = pd.read_table(path, names=["Barcode", "Sample"])
    cells["cell"] = cells["Barcode"] + "." + cells["Sample"]
    return cells


def load_blast_summary(path: Path, level: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    taxon_col = level if level in df.columns else None
    if taxon_col is None:
        for candidate in ("Annotation", "Blast Annotation"):
            if candidate in df.columns:
                taxon_col = candidate
                break
    if taxon_col is None:
        raise SystemExit(
            f"ERROR: no taxon column in {path}; columns={list(df.columns)}"
        )
    out = df.rename(
        columns={
            taxon_col: "taxon",
            "percent": "percent",
            "length": "length",
        }
    )
    needed = {"cell", "taxon", "percent"}
    missing = needed - set(out.columns)
    if missing:
        raise SystemExit(f"ERROR: {path} missing columns {sorted(missing)}")
    out["taxon"] = out["taxon"].astype(str)
    out["percent"] = pd.to_numeric(out["percent"], errors="coerce").fillna(0.0)
    if "length" in out.columns:
        out["length"] = pd.to_numeric(out["length"], errors="coerce").fillna(0.0)
    else:
        out["length"] = 0.0
    out[["Barcode", "Sample"]] = out["cell"].str.split(".", n=1, expand=True)
    return out[["cell", "Barcode", "Sample", "taxon", "percent", "length"]]


def classify_cells(
    blast: pd.DataFrame, candidates: pd.DataFrame, pure_threshold: float
) -> pd.DataFrame:
    """One row per candidate cell: status, top taxon, n taxa, top percent."""
    rows: list[dict[str, object]] = []
    by_cell = {
        cell: group.sort_values("percent", ascending=False)
        for cell, group in blast.groupby("cell", sort=False)
    }
    for row in candidates.itertuples(index=False):
        cell = row.cell
        sample = row.Sample
        group = by_cell.get(cell)
        if group is None or group.empty:
            rows.append(
                {
                    "cell": cell,
                    "Sample": sample,
                    "status": "unannotated",
                    "n_taxa": 0,
                    "top_taxon": "",
                    "top_percent": 0.0,
                    "second_taxon": "",
                    "second_percent": 0.0,
                }
            )
            continue
        top = group.iloc[0]
        second = group.iloc[1] if len(group) > 1 else None
        top_pct = float(top["percent"])
        status = "pure" if top_pct >= pure_threshold else "mixed"
        rows.append(
            {
                "cell": cell,
                "Sample": sample,
                "status": status,
                "n_taxa": int(len(group)),
                "top_taxon": str(top["taxon"]),
                "top_percent": top_pct,
                "second_taxon": "" if second is None else str(second["taxon"]),
                "second_percent": 0.0 if second is None else float(second["percent"]),
            }
        )
    return pd.DataFrame(rows)


def mixed_pair_counts(
    blast: pd.DataFrame, mixed_cells: set[str], min_percent: float
) -> Counter[tuple[str, str]]:
    counts: Counter[tuple[str, str]] = Counter()
    for cell, group in blast.groupby("cell", sort=False):
        if cell not in mixed_cells:
            continue
        taxa = sorted(
            {
                str(taxon)
                for taxon, pct in zip(group["taxon"], group["percent"])
                if float(pct) >= min_percent
            }
        )
        if len(taxa) < 2:
            continue
        for a, b in combinations(taxa, 2):
            counts[(a, b)] += 1
    return counts


def partner_counts(
    blast: pd.DataFrame, mixed_cells: set[str], min_percent: float
) -> tuple[Counter[str], Counter[tuple[str, str]]]:
    """How often a taxon appears in a mixed cell, and as whose partner."""
    in_mixed: Counter[str] = Counter()
    as_partner_of: Counter[tuple[str, str]] = Counter()
    for cell, group in blast.groupby("cell", sort=False):
        if cell not in mixed_cells:
            continue
        ordered = group.sort_values("percent", ascending=False)
        taxa = [
            (str(r.taxon), float(r.percent))
            for r in ordered.itertuples(index=False)
            if float(r.percent) >= min_percent
        ]
        if len(taxa) < 2:
            continue
        for taxon, _pct in taxa:
            in_mixed[taxon] += 1
        dominant = taxa[0][0]
        for taxon, _pct in taxa[1:]:
            as_partner_of[(dominant, taxon)] += 1
    return in_mixed, as_partner_of


def taxon_mixed_share(cells: pd.DataFrame) -> pd.DataFrame:
    """Relative share of mixed cells by top taxon, vs overall abundance.

    Abundant taxa (e.g. Leptocylindrus) naturally appear in more mixed cells;
    enrichment_vs_abundance > 1 means over-represented among mixed tops.
    """
    annotated = cells[cells["status"].isin(["pure", "mixed"])].copy()
    if annotated.empty:
        return pd.DataFrame(
            columns=[
                "taxon",
                "annotated_cells",
                "pure_cells",
                "mixed_as_top",
                "share_of_annotated",
                "share_of_mixed",
                "fraction_mixed_when_top",
                "enrichment_vs_abundance",
            ]
        )
    n_annot = len(annotated)
    n_mixed = int((annotated["status"] == "mixed").sum())
    rows: list[dict[str, object]] = []
    for taxon, group in annotated.groupby("top_taxon", sort=False):
        if not taxon:
            continue
        n = len(group)
        n_pure = int((group["status"] == "pure").sum())
        n_mix = int((group["status"] == "mixed").sum())
        share_annot = n / n_annot if n_annot else 0.0
        share_mixed = n_mix / n_mixed if n_mixed else 0.0
        frac_mixed = n_mix / n if n else 0.0
        enrichment = (share_mixed / share_annot) if share_annot > 0 else 0.0
        rows.append(
            {
                "taxon": str(taxon),
                "annotated_cells": n,
                "pure_cells": n_pure,
                "mixed_as_top": n_mix,
                "share_of_annotated": share_annot,
                "share_of_mixed": share_mixed,
                "fraction_mixed_when_top": frac_mixed,
                "enrichment_vs_abundance": enrichment,
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values(
        ["share_of_mixed", "annotated_cells"], ascending=[False, False]
    ).reset_index(drop=True)


def grouped_share_chart(
    path: Path,
    labels: list[str],
    share_annotated: list[float],
    share_mixed: list[float],
    title: str,
) -> None:
    """Side-by-side horizontal bars: abundance vs share of mixed cells."""
    if not labels:
        path.write_text(
            f'<svg xmlns="http://www.w3.org/2000/svg" width="400" height="80">'
            f'<text x="20" y="40">{_svg_escape(title)}: no data</text></svg>\n'
        )
        return
    n = len(labels)
    width = 960
    height = max(280, 50 + n * 36)
    margin_left, margin_right, margin_top, margin_bottom = 220, 160, 50, 50
    plot_w = width - margin_left - margin_right
    plot_h = height - margin_top - margin_bottom
    vmax = max(share_annotated + share_mixed + [0.01])
    gap = plot_h / max(n, 1)
    bar_h = gap * 0.32
    colors = {"annotated": "#1f4e5f", "mixed": "#c45c26"}
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
        '<rect width="100%" height="100%" fill="#fbfaf6"/>',
        f'<text x="{width // 2}" y="28" text-anchor="middle" '
        f'font-family="Helvetica,Arial,sans-serif" font-size="16" font-weight="700">'
        f"{_svg_escape(title)}</text>",
    ]
    for i, label in enumerate(labels):
        y0 = margin_top + i * gap + gap * 0.15
        for j, (key, value) in enumerate(
            (("annotated", share_annotated[i]), ("mixed", share_mixed[i]))
        ):
            y = y0 + j * (bar_h + 2)
            w = (value / vmax) * plot_w
            parts.append(
                f'<rect x="{margin_left}" y="{y:.1f}" width="{w:.1f}" '
                f'height="{bar_h:.1f}" fill="{colors[key]}"/>'
            )
            parts.append(
                f'<text x="{margin_left + w + 6:.1f}" y="{y + bar_h * 0.78:.1f}" '
                f'font-family="Helvetica,Arial,sans-serif" font-size="10">'
                f"{100 * value:.1f}%</text>"
            )
        parts.append(
            f'<text x="{margin_left - 8}" y="{y0 + bar_h + 4:.1f}" text-anchor="end" '
            f'font-family="Helvetica,Arial,sans-serif" font-size="12">'
            f"{_svg_escape(label)}</text>"
        )
    lx = width - 150
    for j, (key, color) in enumerate(colors.items()):
        label = "share of annotated" if key == "annotated" else "share of mixed"
        parts.append(
            f'<rect x="{lx}" y="{58 + j * 22}" width="14" height="14" fill="{color}"/>'
        )
        parts.append(
            f'<text x="{lx + 20}" y="{69 + j * 22}" '
            f'font-family="Helvetica,Arial,sans-serif" font-size="12">'
            f"{label}</text>"
        )
    parts.append(
        f'<text x="{margin_left + plot_w / 2:.1f}" y="{height - 14}" '
        f'text-anchor="middle" font-family="Helvetica,Arial,sans-serif" '
        f'font-size="12">Relative share</text>'
    )
    parts.append("</svg>\n")
    path.write_text("\n".join(parts))


def write_tsv(path: Path, df: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


# --- SVG helpers (no matplotlib dependency) ---------------------------------

PALETTE = [
    "#1f4e5f",
    "#2e8b57",
    "#c45c26",
    "#3d5a80",
    "#9a031e",
    "#5c4d7a",
    "#0077b6",
    "#6a994e",
    "#bc4749",
    "#588157",
    "#b08968",
    "#4a4e69",
]


def _svg_escape(text: str) -> str:
    return html.escape(str(text), quote=True)


def bar_chart(
    path: Path,
    labels: list[str],
    values: list[float],
    title: str,
    xlabel: str,
    ylabel: str,
    horizontal: bool = True,
    width: int = 900,
    height: int | None = None,
) -> None:
    if not labels:
        path.write_text(
            f'<svg xmlns="http://www.w3.org/2000/svg" width="400" height="80">'
            f'<text x="20" y="40">{_svg_escape(title)}: no data</text></svg>\n'
        )
        return
    n = len(labels)
    if height is None:
        height = max(280, 40 + n * 28) if horizontal else 420
    margin_left = 220 if horizontal else 70
    margin_right = 40
    margin_top = 50
    margin_bottom = 70 if not horizontal else 40
    plot_w = width - margin_left - margin_right
    plot_h = height - margin_top - margin_bottom
    vmax = max(values) if values else 1.0
    vmax = vmax if vmax > 0 else 1.0

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
        f'<rect width="100%" height="100%" fill="#fbfaf6"/>',
        f'<text x="{width // 2}" y="28" text-anchor="middle" '
        f'font-family="Helvetica,Arial,sans-serif" font-size="16" font-weight="700">'
        f"{_svg_escape(title)}</text>",
    ]
    if horizontal:
        bar_h = plot_h / max(n, 1) * 0.7
        gap = plot_h / max(n, 1)
        for i, (label, value) in enumerate(zip(labels, values)):
            y = margin_top + i * gap + (gap - bar_h) / 2
            w = (value / vmax) * plot_w
            color = PALETTE[i % len(PALETTE)]
            parts.append(
                f'<rect x="{margin_left}" y="{y:.1f}" width="{w:.1f}" '
                f'height="{bar_h:.1f}" fill="{color}"/>'
            )
            parts.append(
                f'<text x="{margin_left - 8}" y="{y + bar_h * 0.72:.1f}" '
                f'text-anchor="end" font-family="Helvetica,Arial,sans-serif" '
                f'font-size="12">{_svg_escape(label)}</text>'
            )
            parts.append(
                f'<text x="{margin_left + w + 6:.1f}" y="{y + bar_h * 0.72:.1f}" '
                f'font-family="Helvetica,Arial,sans-serif" font-size="11">'
                f"{value:g}</text>"
            )
        parts.append(
            f'<text x="{margin_left + plot_w / 2:.1f}" y="{height - 12}" '
            f'text-anchor="middle" font-family="Helvetica,Arial,sans-serif" '
            f'font-size="12">{_svg_escape(xlabel)}</text>'
        )
    else:
        bar_w = plot_w / max(n, 1) * 0.7
        gap = plot_w / max(n, 1)
        for i, (label, value) in enumerate(zip(labels, values)):
            x = margin_left + i * gap + (gap - bar_w) / 2
            h = (value / vmax) * plot_h
            y = margin_top + plot_h - h
            color = PALETTE[i % len(PALETTE)]
            parts.append(
                f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w:.1f}" '
                f'height="{h:.1f}" fill="{color}"/>'
            )
            parts.append(
                f'<text x="{x + bar_w / 2:.1f}" y="{height - 40}" '
                f'text-anchor="end" transform="rotate(-40 {x + bar_w / 2:.1f} {height - 40})" '
                f'font-family="Helvetica,Arial,sans-serif" font-size="11">'
                f"{_svg_escape(label)}</text>"
            )
            parts.append(
                f'<text x="{x + bar_w / 2:.1f}" y="{y - 4:.1f}" text-anchor="middle" '
                f'font-family="Helvetica,Arial,sans-serif" font-size="10">'
                f"{value:g}</text>"
            )
        parts.append(
            f'<text x="18" y="{margin_top + plot_h / 2}" '
            f'transform="rotate(-90 18 {margin_top + plot_h / 2})" '
            f'text-anchor="middle" font-family="Helvetica,Arial,sans-serif" '
            f'font-size="12">{_svg_escape(ylabel)}</text>'
        )
    parts.append("</svg>\n")
    path.write_text("\n".join(parts))


def stacked_status_chart(
    path: Path,
    samples: list[str],
    pure: list[int],
    mixed: list[int],
    unannotated: list[int],
    title: str,
) -> None:
    width = max(640, 120 + len(samples) * 90)
    height = 420
    margin_left, margin_right, margin_top, margin_bottom = 70, 40, 50, 80
    plot_w = width - margin_left - margin_right
    plot_h = height - margin_top - margin_bottom
    totals = [p + m + u for p, m, u in zip(pure, mixed, unannotated)]
    vmax = max(totals) if totals else 1
    gap = plot_w / max(len(samples), 1)
    bar_w = gap * 0.55
    colors = {"pure": "#2e8b57", "mixed": "#c45c26", "unannotated": "#9aa0a6"}
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
        '<rect width="100%" height="100%" fill="#fbfaf6"/>',
        f'<text x="{width // 2}" y="28" text-anchor="middle" '
        f'font-family="Helvetica,Arial,sans-serif" font-size="16" font-weight="700">'
        f"{_svg_escape(title)}</text>",
    ]
    for i, sample in enumerate(samples):
        x = margin_left + i * gap + (gap - bar_w) / 2
        y = margin_top + plot_h
        for key, value in (
            ("unannotated", unannotated[i]),
            ("mixed", mixed[i]),
            ("pure", pure[i]),
        ):
            h = (value / vmax) * plot_h if vmax else 0
            y -= h
            parts.append(
                f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w:.1f}" '
                f'height="{h:.1f}" fill="{colors[key]}"/>'
            )
        parts.append(
            f'<text x="{x + bar_w / 2:.1f}" y="{height - 45}" text-anchor="end" '
            f'transform="rotate(-35 {x + bar_w / 2:.1f} {height - 45})" '
            f'font-family="Helvetica,Arial,sans-serif" font-size="11">'
            f"{_svg_escape(sample)}</text>"
        )
    # legend
    lx = width - 180
    for j, (key, color) in enumerate(colors.items()):
        parts.append(
            f'<rect x="{lx}" y="{60 + j * 22}" width="14" height="14" fill="{color}"/>'
        )
        parts.append(
            f'<text x="{lx + 20}" y="{71 + j * 22}" '
            f'font-family="Helvetica,Arial,sans-serif" font-size="12">'
            f"{key}</text>"
        )
    parts.append("</svg>\n")
    path.write_text("\n".join(parts))


def heatmap(
    path: Path,
    labels: list[str],
    matrix: list[list[float]],
    title: str,
) -> None:
    n = len(labels)
    if n == 0:
        path.write_text(
            f'<svg xmlns="http://www.w3.org/2000/svg" width="400" height="80">'
            f'<text x="20" y="40">{_svg_escape(title)}: no data</text></svg>\n'
        )
        return
    cell = 28
    margin_left = 160
    margin_top = 140
    width = margin_left + n * cell + 40
    height = margin_top + n * cell + 40
    flat = [v for row in matrix for v in row]
    vmax = max(flat) if flat else 1.0
    vmax = vmax if vmax > 0 else 1.0

    def color(value: float) -> str:
        if value <= 0:
            return "#f0efe8"
        t = math.sqrt(value / vmax)
        r = int(251 - t * 180)
        g = int(250 - t * 200)
        b = int(246 - t * 120)
        return f"rgb({r},{g},{b})"

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
        '<rect width="100%" height="100%" fill="#fbfaf6"/>',
        f'<text x="{width // 2}" y="28" text-anchor="middle" '
        f'font-family="Helvetica,Arial,sans-serif" font-size="16" font-weight="700">'
        f"{_svg_escape(title)}</text>",
    ]
    for i, label in enumerate(labels):
        parts.append(
            f'<text x="{margin_left - 8}" y="{margin_top + i * cell + cell * 0.68:.1f}" '
            f'text-anchor="end" font-family="Helvetica,Arial,sans-serif" font-size="10">'
            f"{_svg_escape(label)}</text>"
        )
        parts.append(
            f'<text x="{margin_left + i * cell + cell / 2:.1f}" y="{margin_top - 8}" '
            f'text-anchor="start" transform="rotate(-60 {margin_left + i * cell + cell / 2:.1f} {margin_top - 8})" '
            f'font-family="Helvetica,Arial,sans-serif" font-size="10">'
            f"{_svg_escape(label)}</text>"
        )
    for i in range(n):
        for j in range(n):
            value = matrix[i][j]
            x = margin_left + j * cell
            y = margin_top + i * cell
            parts.append(
                f'<rect x="{x}" y="{y}" width="{cell - 1}" height="{cell - 1}" '
                f'fill="{color(value)}"/>'
            )
            if value > 0:
                parts.append(
                    f'<text x="{x + cell / 2:.1f}" y="{y + cell * 0.65:.1f}" '
                    f'text-anchor="middle" font-family="Helvetica,Arial,sans-serif" '
                    f'font-size="9">{value:g}</text>'
                )
    parts.append("</svg>\n")
    path.write_text("\n".join(parts))


def write_html(
    path: Path,
    title: str,
    summary_lines: list[str],
    figures: list[tuple[str, Path]],
) -> None:
    body = ["<html><head><meta charset='utf-8'>", f"<title>{html.escape(title)}</title>"]
    body.append(
        "<style>body{font-family:Helvetica,Arial,sans-serif;margin:24px;background:#fbfaf6;"
        "color:#1b1b1b} pre{background:#fff;padding:12px;border:1px solid #ddd}"
        "img,object{max-width:100%;margin:12px 0}</style></head><body>"
    )
    body.append(f"<h1>{html.escape(title)}</h1>")
    body.append("<pre>" + html.escape("\n".join(summary_lines)) + "</pre>")
    for caption, fig in figures:
        body.append(f"<h2>{html.escape(caption)}</h2>")
        body.append(f"<object type='image/svg+xml' data='{html.escape(fig.name)}'></object>")
    body.append("</body></html>\n")
    path.write_text("\n".join(body))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--blast-summary", required=True)
    parser.add_argument("--summarized", required=True)
    parser.add_argument("--candidates", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--level", required=True)
    parser.add_argument("--chimeric-percent", type=float, default=80.0)
    parser.add_argument(
        "--mixed-min-percent",
        type=float,
        default=5.0,
        help="Minimum taxon percent (0-100) to count in mixed-pair analysis.",
    )
    parser.add_argument("--top-taxa", type=int, default=20)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    pure_threshold = args.chimeric_percent / 100.0
    mixed_min = args.mixed_min_percent / 100.0

    candidates = load_candidates(Path(args.candidates))
    blast = load_blast_summary(Path(args.blast_summary), args.level)
    cells = classify_cells(blast, candidates, pure_threshold)
    mixed_cells = set(cells.loc[cells["status"] == "mixed", "cell"])
    summarized = pd.read_csv(Path(args.summarized))
    summarized_cells = (
        set(summarized["cell"].astype(str)) if "cell" in summarized.columns else set()
    )

    # Overview by sample
    samples = list(dict.fromkeys(candidates["Sample"].tolist()))
    overview_rows = []
    for sample in samples + ["All"]:
        subset = cells if sample == "All" else cells[cells["Sample"] == sample]
        n = len(subset)
        n_pure = int((subset["status"] == "pure").sum())
        n_mixed = int((subset["status"] == "mixed").sum())
        n_un = int((subset["status"] == "unannotated").sum())
        n_annot = n_pure + n_mixed
        overview_rows.append(
            {
                "sample": sample,
                "cells": n,
                "annotated": n_annot,
                "pure": n_pure,
                "mixed": n_mixed,
                "unannotated": n_un,
                "fraction_annotated": (n_annot / n) if n else 0.0,
                "fraction_pure_of_all": (n_pure / n) if n else 0.0,
                "fraction_pure_of_annotated": (n_pure / n_annot) if n_annot else 0.0,
                "fraction_mixed_of_annotated": (n_mixed / n_annot) if n_annot else 0.0,
            }
        )
    overview = pd.DataFrame(overview_rows)
    write_tsv(outdir / "host_18s_qc_overview.tsv", overview)

    # Common taxa among pure cells
    pure_cells = cells[cells["status"] == "pure"]
    taxa_counts = pure_cells["top_taxon"].value_counts().reset_index()
    taxa_counts.columns = ["taxon", "pure_cells"]
    taxa_counts["fraction_of_pure"] = taxa_counts["pure_cells"] / max(len(pure_cells), 1)
    # also by sample
    by_sample = (
        pure_cells.groupby(["Sample", "top_taxon"])
        .size()
        .reset_index(name="pure_cells")
        .sort_values(["Sample", "pure_cells"], ascending=[True, False])
    )
    write_tsv(outdir / "host_18s_qc_taxa.tsv", taxa_counts)
    write_tsv(outdir / "host_18s_qc_taxa_by_sample.tsv", by_sample)

    # Mixed pairs and partners
    pair_counts = mixed_pair_counts(blast, mixed_cells, mixed_min)
    pairs_df = pd.DataFrame(
        [
            {"taxon_a": a, "taxon_b": b, "mixed_cells": n}
            for (a, b), n in pair_counts.most_common()
        ]
    )
    if pairs_df.empty:
        pairs_df = pd.DataFrame(columns=["taxon_a", "taxon_b", "mixed_cells"])
    write_tsv(outdir / "host_18s_qc_mixed_pairs.tsv", pairs_df)

    in_mixed, as_partner = partner_counts(blast, mixed_cells, mixed_min)
    partners_df = pd.DataFrame(
        [
            {
                "dominant_taxon": dom,
                "partner_taxon": partner,
                "mixed_cells": n,
            }
            for (dom, partner), n in as_partner.most_common()
        ]
    )
    if partners_df.empty:
        partners_df = pd.DataFrame(
            columns=["dominant_taxon", "partner_taxon", "mixed_cells"]
        )
    write_tsv(outdir / "host_18s_qc_mixed_partners.tsv", partners_df)

    n_mixed_total = len(mixed_cells)
    mixed_freq = pd.DataFrame(
        [{"taxon": t, "mixed_cells": n} for t, n in in_mixed.most_common()]
    )
    if mixed_freq.empty:
        mixed_freq = pd.DataFrame(
            columns=["taxon", "mixed_cells", "share_of_mixed"]
        )
    else:
        mixed_freq["share_of_mixed"] = mixed_freq["mixed_cells"] / max(
            n_mixed_total, 1
        )
    write_tsv(outdir / "host_18s_qc_taxa_in_mixed.tsv", mixed_freq)

    # Relative share of mixed cells by top taxon vs overall abundance
    share_df = taxon_mixed_share(cells)
    write_tsv(outdir / "host_18s_qc_mixed_share.tsv", share_df)

    write_tsv(outdir / "host_18s_qc_cells.tsv", cells)

    # Figures
    top_n = max(1, args.top_taxa)
    top_taxa = taxa_counts.head(top_n)
    fig_taxa = outdir / "host_18s_qc_taxa_bar.svg"
    bar_chart(
        fig_taxa,
        top_taxa["taxon"].tolist(),
        top_taxa["pure_cells"].astype(float).tolist(),
        title=f"Top {args.level} among pure cells",
        xlabel="Pure cells",
        ylabel=args.level,
        horizontal=True,
    )

    fig_status = outdir / "host_18s_qc_status_by_sample.svg"
    stacked_status_chart(
        fig_status,
        samples,
        [int(overview.loc[overview["sample"] == s, "pure"].iloc[0]) for s in samples],
        [int(overview.loc[overview["sample"] == s, "mixed"].iloc[0]) for s in samples],
        [
            int(overview.loc[overview["sample"] == s, "unannotated"].iloc[0])
            for s in samples
        ],
        title="Host 18S annotation status by sample",
    )

    top_pairs = pairs_df.head(top_n)
    fig_pairs = outdir / "host_18s_qc_mixed_pairs.svg"
    pair_labels = [
        f"{a} + {b}" for a, b in zip(top_pairs.get("taxon_a", []), top_pairs.get("taxon_b", []))
    ]
    bar_chart(
        fig_pairs,
        pair_labels,
        top_pairs["mixed_cells"].astype(float).tolist() if not top_pairs.empty else [],
        title=f"Most common mixed {args.level} pairs (≥{args.mixed_min_percent:g}% each)",
        xlabel="Mixed cells",
        ylabel="Pair",
        horizontal=True,
    )

    # Co-occurrence heatmap for taxa that appear most often in mixed cells
    heat_taxa = mixed_freq["taxon"].head(min(15, top_n)).tolist()
    if not heat_taxa and not taxa_counts.empty:
        heat_taxa = taxa_counts["taxon"].head(min(15, top_n)).tolist()
    index = {t: i for i, t in enumerate(heat_taxa)}
    matrix = [[0.0 for _ in heat_taxa] for _ in heat_taxa]
    for (a, b), n in pair_counts.items():
        if a in index and b in index:
            matrix[index[a]][index[b]] = float(n)
            matrix[index[b]][index[a]] = float(n)
    fig_heat = outdir / "host_18s_qc_cooccurrence.svg"
    heatmap(
        fig_heat,
        heat_taxa,
        matrix,
        title=f"Co-occurrence of {args.level} in mixed cells",
    )

    fig_partners = outdir / "host_18s_qc_mixed_taxa.svg"
    bar_chart(
        fig_partners,
        mixed_freq["taxon"].head(top_n).tolist(),
        mixed_freq["mixed_cells"].head(top_n).astype(float).tolist()
        if not mixed_freq.empty
        else [],
        title=f"{args.level} most often present in mixed cells",
        xlabel="Mixed cells containing taxon",
        ylabel=args.level,
        horizontal=True,
    )

    top_share = share_df.head(top_n)
    fig_share = outdir / "host_18s_qc_mixed_share.svg"
    grouped_share_chart(
        fig_share,
        top_share["taxon"].tolist() if not top_share.empty else [],
        top_share["share_of_annotated"].astype(float).tolist()
        if not top_share.empty
        else [],
        top_share["share_of_mixed"].astype(float).tolist()
        if not top_share.empty
        else [],
        title=f"Relative share of mixed cells vs abundance ({args.level})",
    )

    # Sort by mixed-when-top so the figure highlights taxa that are rarely pure
    mixed_when_top = (
        share_df.sort_values("fraction_mixed_when_top", ascending=False).head(top_n)
        if not share_df.empty
        else share_df
    )
    fig_mixed_when_top = outdir / "host_18s_qc_mixed_when_top.svg"
    bar_chart(
        fig_mixed_when_top,
        mixed_when_top["taxon"].tolist() if not mixed_when_top.empty else [],
        [
            100.0 * float(v)
            for v in (
                mixed_when_top["fraction_mixed_when_top"].tolist()
                if not mixed_when_top.empty
                else []
            )
        ],
        title=f"Mixed when top ({args.level}): % of cells with this top taxon that are mixed",
        xlabel="% mixed when top",
        ylabel=args.level,
        horizontal=True,
    )

    # Text report
    all_row = overview[overview["sample"] == "All"].iloc[0]
    lines = [
        f"Host 18S QC at taxonomic level: {args.level}",
        f"Pure-cell threshold: top taxon ≥ {args.chimeric_percent:g}% of 18S alignment length",
        f"Mixed-pair floor: each taxon ≥ {args.mixed_min_percent:g}%",
        "",
        f"Cells: {int(all_row['cells'])}",
        (
            f"Annotated: {int(all_row['annotated'])} "
            f"({100 * float(all_row['fraction_annotated']):.1f}% of cells)"
        ),
        (
            f"Pure: {int(all_row['pure'])} "
            f"({100 * float(all_row['fraction_pure_of_all']):.1f}% of cells, "
            f"{100 * float(all_row['fraction_pure_of_annotated']):.1f}% of annotated)"
        ),
        (
            f"Mixed: {int(all_row['mixed'])} "
            f"({100 * float(all_row['fraction_mixed_of_annotated']):.1f}% of annotated)"
        ),
        f"Unannotated: {int(all_row['unannotated'])}",
        (
            f"Summarized pure table: {len(summarized_cells)} cells "
            f"(QC pure={int(all_row['pure'])})"
        ),
        "",
        "Most common taxa among pure cells:",
    ]
    for row in taxa_counts.head(15).itertuples(index=False):
        lines.append(
            f"  {row.taxon}: {int(row.pure_cells)} "
            f"({100 * float(row.fraction_of_pure):.1f}% of pure)"
        )
    lines.append("")
    lines.append(
        "Relative share of mixed cells by top taxon "
        "(abundant taxa dominate mixed cells by chance; "
        "enrichment>1 = over-represented vs abundance):"
    )
    if share_df.empty:
        lines.append("  (none)")
    else:
        for row in share_df.head(15).itertuples(index=False):
            lines.append(
                f"  {row.taxon}: {100 * float(row.share_of_mixed):.1f}% of mixed "
                f"(vs {100 * float(row.share_of_annotated):.1f}% of annotated; "
                f"mixed_when_top={100 * float(row.fraction_mixed_when_top):.1f}%; "
                f"enrichment={float(row.enrichment_vs_abundance):.2f})"
            )
    lines.append("")
    lines.append("Most common mixed pairs:")
    if pairs_df.empty:
        lines.append("  (none)")
    else:
        for row in pairs_df.head(15).itertuples(index=False):
            lines.append(f"  {row.taxon_a} + {row.taxon_b}: {int(row.mixed_cells)} cells")
    lines.append("")
    lines.append("Taxa most often present in mixed cells:")
    if mixed_freq.empty:
        lines.append("  (none)")
    else:
        for row in mixed_freq.head(15).itertuples(index=False):
            share = getattr(row, "share_of_mixed", 0.0)
            lines.append(
                f"  {row.taxon}: {int(row.mixed_cells)} mixed cells "
                f"({100 * float(share):.1f}% of mixed)"
            )
    lines.append("")
    lines.append("By sample:")
    for row in overview.itertuples(index=False):
        if row.sample == "All":
            continue
        lines.append(
            f"  {row.sample}: pure={int(row.pure)} mixed={int(row.mixed)} "
            f"unannotated={int(row.unannotated)} "
            f"(pure/annotated={100 * float(row.fraction_pure_of_annotated):.1f}%)"
        )

    report_path = outdir / "host_18s_qc.txt"
    report_path.write_text("\n".join(lines) + "\n")

    write_html(
        outdir / "host_18s_qc.html",
        f"Host 18S QC ({args.level})",
        lines,
        [
            ("Annotation status by sample", fig_status),
            (f"Top {args.level} among pure cells", fig_taxa),
            ("Relative share of mixed vs abundance", fig_share),
            ("Mixed when top", fig_mixed_when_top),
            ("Common mixed pairs", fig_pairs),
            ("Taxa in mixed cells", fig_partners),
            ("Co-occurrence in mixed cells", fig_heat),
        ],
    )
    print(f"[host_18s_qc] wrote {report_path}", flush=True)
    print(f"[host_18s_qc] wrote {outdir / 'host_18s_qc.html'}", flush=True)


if __name__ == "__main__":
    main()
