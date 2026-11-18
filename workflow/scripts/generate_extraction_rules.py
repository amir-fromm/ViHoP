"""Generate static Snakemake rules for per-sample read extraction."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

from common import _write_tsv_if_changed, cellranger_bam, load_candidate_cells, resolve_bam_root


def write_extraction_rules(
    output_path: str | Path,
    samples: list[str],
    config: dict,
    pending_extract: str | Path,
) -> None:
    """Write one extract_reads rule per sample.

    The rule input is results/tmp_extract_<sample>.tsv, not the full candidate
    list. Samples with no pending cells get a local rule that only stamps
    results/extraction/<sample>.done, so the BAM is not scanned again.
    """
    resolve_bam_root(config)  # validate / document default
    outdir = config["outputs"]["cells_dir"]
    threads = config["params"]["extraction"]["threads"]
    mem_mb = config["resources"]["extraction"]["mem_mb"]
    disk_mb = config["resources"]["extraction"].get("disk_mb", 4000)
    runtime = config["resources"]["extraction"].get("runtime", 720)
    lsf_queue = config["resources"]["extraction"].get(
        "lsf_queue", config.get("_lsf_queue_short", "short")
    )
    extract_script = str(
        Path(__file__).resolve().parent / "extract_reads_from_bam.py"
    )
    pending_path = Path(pending_extract)
    grouped: dict[str, list[tuple[str, str]]] = defaultdict(list)
    if pending_path.is_file():
        for barcode, sample in load_candidate_cells(pending_path):
            grouped[sample].append((barcode, sample))

    blocks: list[str] = []
    local_rules: list[str] = []
    for sample in samples:
        # One done-flag per sample. Listing every cell FASTQ here makes Snakemake
        # stat tens of thousands of files while building the DAG.
        done = f"results/extraction/{sample}.done"
        cells = f"results/tmp_extract_{sample}.tsv"
        _write_tsv_if_changed(Path(cells), grouped.get(sample, []))
        if not grouped.get(sample):
            local_rules.append(f"extract_reads_{sample}")
            blocks.append(
                f"""rule extract_reads_{sample}:
    input:
        cells="{cells}",
    output:
        done="{done}",
    shell:
        \"\"\"
        mkdir -p "$(dirname {{output.done}})"
        echo "[extract] {sample}: every candidate cell already has reads.fq.gz" > {{output.done}}
        \"\"\"
"""
            )
            continue
        bam = cellranger_bam(sample, config)
        blocks.append(
            f"""rule extract_reads_{sample}:
    input:
        bam="{bam}",
        cells="{cells}",
        script="{extract_script}",
    output:
        done="{done}",
    params:
        outdir="{outdir}",
    threads: {threads}
    resources:
        mem_mb={mem_mb},
        disk_mb={disk_mb},
        runtime={runtime},
        lsf_queue="{lsf_queue}",
    shell:
        \"\"\"
        python {{input.script}} \\
            --bam {{input.bam}} \\
            --cells {{input.cells}} \\
            --outdir {{params.outdir}} \\
            --threads {{threads}}
        mkdir -p "$(dirname {{output.done}})"
        touch {{output.done}}
        \"\"\"
"""
        )

    text = ""
    if local_rules:
        text += "localrules:\n    " + ",\n    ".join(local_rules) + "\n\n"
    text += "\n".join(blocks)
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.is_file() or path.read_text() != text:
        path.write_text(text)
