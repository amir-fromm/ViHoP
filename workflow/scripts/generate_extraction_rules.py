"""Generate static Snakemake rules for per-sample read extraction."""

from __future__ import annotations

from pathlib import Path

from common import cell_outputs_for_sample


def write_extraction_rules(
    output_path: str | Path,
    candidate_cells: str | Path,
    samples: list[str],
    config: dict,
) -> None:
    """Write one traditional extract_reads rule per sample."""
    bam_root = config["references"]["bam_root"]
    outdir = config["outputs"]["cells_dir"]
    threads = config["params"]["extraction"]["threads"]
    mem_mb = config["resources"]["extraction"]["mem_mb"]
    disk_mb = config["resources"]["extraction"].get("disk_mb", 4000)
    runtime = config["resources"]["extraction"].get("runtime", 720)
    lsf_queue = config["resources"]["extraction"].get("lsf_queue", "short")

    blocks: list[str] = []
    for sample in samples:
        outputs = cell_outputs_for_sample(candidate_cells, sample)
        output_lines = ",\n        ".join(f'"{path}"' for path in outputs)
        blocks.append(
            f"""rule extract_reads_{sample}:
    input:
        bam="{bam_root}/{sample}/outs/possorted_genome_bam.bam",
        cells="results/candidate_cells_{sample}.tsv",
    output:
        {output_lines},
    params:
        outdir="{outdir}",
    threads: {threads},
    resources:
        mem_mb={mem_mb},
        disk_mb={disk_mb},
        runtime={runtime},
        lsf_queue="{lsf_queue}",
    shell:
        \"\"\"
        python workflow/scripts/extract_reads_from_bam.py \\
            --bam {{input.bam}} \\
            --cells {{input.cells}} \\
            --outdir {{params.outdir}} \\
            --threads {{threads}}
        \"\"\"
"""
        )

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(blocks))
