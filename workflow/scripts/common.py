"""Shared helpers for the ViHoP Snakemake workflow."""

from __future__ import annotations

from pathlib import Path


def load_candidate_cells(path: str | Path) -> list[tuple[str, str]]:
    """Load barcode/sample pairs from a tab-delimited candidate cell list."""
    cells: list[tuple[str, str]] = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            barcode, sample = line.split("\t", 1)
            cells.append((barcode, sample))
    return cells


def cell_id(barcode: str, sample: str) -> str:
    return f"{barcode}.{sample}"


def split_cell_id(cell: str) -> tuple[str, str]:
    barcode, sample = cell.rsplit(".", 1)
    return barcode, sample


def get_all_cell_ids(path: str | Path) -> list[str]:
    return [cell_id(barcode, sample) for barcode, sample in load_candidate_cells(path)]


def get_samples(path: str | Path) -> list[str]:
    return sorted({sample for _, sample in load_candidate_cells(path)})


def get_cells_by_sample(path: str | Path) -> dict[str, list[str]]:
    grouped: dict[str, list[str]] = {}
    for barcode, sample in load_candidate_cells(path):
        grouped.setdefault(sample, []).append(cell_id(barcode, sample))
    return grouped


def cell_reads(barcode: str, sample: str) -> str:
    return f"results/cells/{barcode}.{sample}/reads.fq.gz"


def cell_metadata(barcode: str, sample: str) -> str:
    return f"results/cells/{barcode}.{sample}/metadata.yaml"


def cell_transcripts(barcode: str, sample: str) -> str:
    return f"results/cells/{barcode}.{sample}/transcripts.fasta"


def cell_outputs_for_sample(path: str | Path, sample: str) -> list[str]:
    outputs: list[str] = []
    for barcode, cell_sample in load_candidate_cells(path):
        if cell_sample != sample:
            continue
        outputs.append(cell_reads(barcode, cell_sample))
        outputs.append(cell_metadata(barcode, cell_sample))
    return outputs


def all_cell_extraction_outputs(path: str | Path) -> list[str]:
    outputs: list[str] = []
    for barcode, sample in load_candidate_cells(path):
        outputs.append(cell_reads(barcode, sample))
        outputs.append(cell_metadata(barcode, sample))
    return outputs


def all_cell_transcripts(path: str | Path) -> list[str]:
    return [
        cell_transcripts(barcode, sample)
        for barcode, sample in load_candidate_cells(path)
    ]


def cell_transcripts_edit(barcode: str, sample: str) -> str:
    return f"results/cells/{barcode}.{sample}/transcripts.edit.fasta"


def all_cell_transcripts_edit(path: str | Path) -> list[str]:
    return [
        cell_transcripts_edit(barcode, sample)
        for barcode, sample in load_candidate_cells(path)
    ]


def cell_lastal_18s(barcode: str, sample: str) -> str:
    return f"results/host/lastal/{barcode}.{sample}.last"


def all_cell_lastal_18s(path: str | Path) -> list[str]:
    return [
        cell_lastal_18s(barcode, sample)
        for barcode, sample in load_candidate_cells(path)
    ]
