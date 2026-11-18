"""Shared helpers for the ViHoP Snakemake workflow."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path


def load_sample_table(path: str | Path) -> dict[str, list[str]]:
    """Load sample_id -> list of FASTQ directories from a TSV.

    Column 2 may be a single path or comma-separated paths (multi-run / reseq).
    """
    samples: dict[str, list[str]] = {}
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            sample, fastqs = line.split("\t", 1)
            paths = [p.strip() for p in fastqs.split(",") if p.strip()]
            if not paths:
                raise ValueError(f"No FASTQ paths for sample {sample!r} in {path}")
            samples[sample] = paths
    return samples


def resolve_cellranger_samples(config: dict) -> list[str]:
    """Sample IDs for Cell Ranger: config input.samples, or all rows in sample_table."""
    table = load_sample_table(config["input"]["sample_table"])
    selected = config["input"].get("samples")
    if not selected:
        return sorted(table)
    missing = [s for s in selected if s not in table]
    if missing:
        raise ValueError(
            f"input.samples not found in sample_table: {missing}. "
            f"Known: {sorted(table)}"
        )
    return list(selected)


def sample_fastq_paths(config: dict, sample: str) -> list[str]:
    return load_sample_table(config["input"]["sample_table"])[sample]


def sample_fastqs_csv(config: dict, sample: str) -> str:
    """Comma-separated FASTQ dirs for cellranger --fastqs."""
    return ",".join(sample_fastq_paths(config, sample))


_DEFAULT_BAM_ROOT = "results/cellranger"


def resolve_bam_root(config: dict) -> str:
    """Parent of <sample>/outs/... Default: results/cellranger (under workdir).

    Override via config references.bam_root or ``./vihop --bam-root PATH``.
    """
    root = (config.get("references") or {}).get("bam_root", None)
    if root in (None, "", "null", False, "false"):
        return _DEFAULT_BAM_ROOT
    return str(root).rstrip("/")


def cellranger_sample_outdir(sample: str, config: dict | None = None) -> str:
    """Directory that contains outs/ for a sample (under bam_root)."""
    root = resolve_bam_root(config) if config is not None else _DEFAULT_BAM_ROOT
    return f"{root}/{sample}"


def cellranger_bam(sample: str, config: dict | None = None) -> str:
    return f"{cellranger_sample_outdir(sample, config)}/outs/possorted_genome_bam.bam"


def cellranger_filtered_matrix(sample: str, config: dict | None = None) -> str:
    return f"{cellranger_sample_outdir(sample, config)}/outs/filtered_feature_bc_matrix"


def cellranger_barcodes(sample: str, config: dict | None = None) -> str:
    """Marker file inside the filtered matrix (choose_cells does not need the BAM)."""
    return f"{cellranger_filtered_matrix(sample, config)}/barcodes.tsv.gz"


def all_cellranger_bams(config: dict) -> list[str]:
    return [cellranger_bam(s, config) for s in resolve_cellranger_samples(config)]


def all_cellranger_barcodes(config: dict) -> list[str]:
    return [cellranger_barcodes(s, config) for s in resolve_cellranger_samples(config)]


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


def _write_tsv_if_changed(path: Path, rows: list[tuple[str, str]]) -> None:
    text = "".join(f"{barcode}\t{sample}\n" for barcode, sample in rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.is_file() and path.read_text() == text:
        return
    path.write_text(text)


def _nonempty_file(path: Path) -> bool:
    try:
        return path.is_file() and path.stat().st_size > 0
    except OSError:
        return False


def write_pending_cell_lists(
    candidates: str | Path,
    cells_dir: str | Path,
    extract_path: str | Path = "results/tmp_extract.tsv",
    assemble_path: str | Path = "results/tmp_assemble.tsv",
    stamp_path: str | Path = "results/tmp_pending.ready",
) -> tuple[list[tuple[str, str]], list[tuple[str, str]]]:
    """Write cells that still need extraction or assembly.

    A cell is omitted from tmp_extract when reads.fq.gz is already non-empty.
    It is omitted from tmp_assemble when transcripts.fasta and
    transcripts.edit.fasta are both non-empty. Cells that will be extracted
    are always assembled, because extraction replaces their transcripts.
    Existing TSV bytes are left untouched when the rows are unchanged, so
    Snakemake does not rerun a finished sample just to rewrite the list.
    """
    rows = load_candidate_cells(candidates)
    cells = Path(cells_dir)
    need_extract: list[tuple[str, str]] = []
    need_assemble: list[tuple[str, str]] = []
    for barcode, sample in rows:
        folder = cells / cell_id(barcode, sample)
        has_reads = _nonempty_file(folder / "reads.fq.gz")
        has_transcripts = _nonempty_file(folder / "transcripts.fasta")
        has_edit = _nonempty_file(folder / "transcripts.edit.fasta")
        if not has_reads:
            need_extract.append((barcode, sample))
        if (not has_reads) or (not has_transcripts) or (not has_edit):
            need_assemble.append((barcode, sample))

    extract_file = Path(extract_path)
    _write_tsv_if_changed(extract_file, need_extract)
    _write_tsv_if_changed(Path(assemble_path), need_assemble)

    by_sample: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for barcode, sample in need_extract:
        by_sample[sample].append((barcode, sample))
    samples = sorted({sample for _, sample in rows})
    for sample in samples:
        sample_path = extract_file.parent / f"tmp_extract_{sample}.tsv"
        _write_tsv_if_changed(sample_path, by_sample.get(sample, []))

    stamp = Path(stamp_path)
    stamp.parent.mkdir(parents=True, exist_ok=True)
    stamp.write_text(
        f"extract={len(need_extract)}\nassemble={len(need_assemble)}\n"
    )
    return need_extract, need_assemble


def pending_lists_ready(
    extract_path: str | Path = "results/tmp_extract.tsv",
    assemble_path: str | Path = "results/tmp_assemble.tsv",
    stamp_path: str | Path = "results/tmp_pending.ready",
) -> bool:
    """True when this run already wrote the pending lists.

    Cluster jobs re-parse the Snakefile. They must reuse the lists from the
    controller instead of scanning results/cells again while jobs are writing.
    """
    return all(Path(path).is_file() for path in (extract_path, assemble_path, stamp_path))


def delete_pending_cell_lists(results_dir: str | Path) -> None:
    """Remove tmp_extract, tmp_assemble, and the per-sample extract lists."""
    results = Path(results_dir)
    for name in ("tmp_extract.tsv", "tmp_assemble.tsv", "tmp_pending.ready"):
        path = results / name
        if path.is_file():
            path.unlink()
    if results.is_dir():
        for path in results.glob("tmp_extract_*.tsv"):
            if path.is_file():
                path.unlink()


def write_assembly_batches(
    candidate_cells: str | Path,
    batch_size: int,
    outdir: str | Path = "results/assembly/batches",
) -> list[str]:
    """Write cell lists for batched assembly. Returns batch names (batch_0000, ...)."""
    rows = load_candidate_cells(candidate_cells)
    if batch_size < 1:
        raise ValueError(f"assembly batch_size must be >= 1, got {batch_size}")
    dest = Path(outdir)
    dest.mkdir(parents=True, exist_ok=True)
    names: list[str] = []
    for index in range(0, len(rows), batch_size):
        name = f"batch_{index // batch_size:04d}"
        chunk = rows[index : index + batch_size]
        text = "".join(f"{barcode}\t{sample}\n" for barcode, sample in chunk)
        path = dest / f"{name}.tsv"
        if not path.is_file() or path.read_text() != text:
            path.write_text(text)
        names.append(name)
    return names


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


def apply_lsf_queue_names(config: dict) -> dict:
    """Map lsf_queue roles onto the queue names of this cluster.

    resources.*.lsf_queue is a role, usually ``short`` or ``long``.
    ``lsf.queues.<role>`` is the name the scheduler actually has.
    A value that is not a role is kept as written.
    """
    raw = ((config.get("lsf") or {}).get("queues") or {})
    mapping: dict[str, str] = {}
    for role, name in raw.items():
        text = "" if name is None else str(name).strip()
        if text:
            mapping[str(role)] = text
    for group in (config.get("resources") or {}).values():
        if not isinstance(group, dict) or "lsf_queue" not in group:
            continue
        role = str(group["lsf_queue"]).strip()
        if role in mapping:
            group["lsf_queue"] = mapping[role]
    config["_lsf_queue_short"] = mapping.get("short", "short")
    config["_lsf_queue_long"] = mapping.get("long", "long")
    return config


def load_finished_steps(path: str | Path = "results/finished/steps.yaml") -> dict[str, bool]:
    """Read results/finished/steps.yaml. Missing file means nothing is finished."""
    steps: dict[str, bool] = {}
    file = Path(path)
    if not file.is_file():
        return steps
    for line in file.read_text().splitlines():
        text = line.split("#", 1)[0].strip()
        if not text or ":" not in text:
            continue
        key, value = text.split(":", 1)
        steps[key.strip()] = value.strip().lower() in {"1", "true", "yes", "done"}
    return steps


def mark_step_finished(name: str, path: str | Path = "results/finished/steps.yaml") -> None:
    """Record one finished workflow step without dropping the other entries."""
    file = Path(path)
    steps = load_finished_steps(file)
    steps[name] = True
    file.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"{key}: {'true' if steps[key] else 'false'}" for key in sorted(steps)]
    file.write_text("\n".join(lines) + "\n")


def cell_lastal_18s(barcode: str, sample: str) -> str:
    return f"results/host/lastal/{barcode}.{sample}.last"


def all_cell_lastal_18s(path: str | Path) -> list[str]:
    return [
        cell_lastal_18s(barcode, sample)
        for barcode, sample in load_candidate_cells(path)
    ]
