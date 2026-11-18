# Build results/candidate_cells.tsv from Cell Ranger filtered matrices.
# Default: all barcodes in filtered_feature_bc_matrix.
# Optional: params.choose_cells.min_umi / min_genes (null → no extra filter).

from common import all_cellranger_barcodes, resolve_bam_root, resolve_cellranger_samples

CC_SAMPLES = resolve_cellranger_samples(config)
CC_CFG = config["params"]["choose_cells"]
CC_RES = config["resources"].get("choose_cells", {})
BAM_ROOT = resolve_bam_root(config)
CC_OUT = config["outputs"]["candidate_cells"]
CC_SCRIPT = vihop_script("workflow/scripts/run_choose_cells.py")

# null / false / 0 → omit filter (keep all filtered cells)
_raw_umi = CC_CFG.get("min_umi", None)
_raw_genes = CC_CFG.get("min_genes", CC_CFG.get("min_gene", None))


def _filter_arg(val) -> str:
    if val in (None, "null", False, "false", 0, "0"):
        return ""
    return str(int(val))


CC_MIN_UMI = _filter_arg(_raw_umi)
CC_MIN_GENES = _filter_arg(_raw_genes)
CC_MEM = CC_RES.get("mem_mb", 8000)
CC_DISK = CC_RES.get("disk_mb", 8000)
CC_RUNTIME = CC_RES.get("runtime", 120)
CC_QUEUE = CC_RES.get("lsf_queue", config.get("_lsf_queue_short", "short"))


rule choose_cells:
    input:
        # Filtered matrix only (BAM not required). barcodes.tsv.gz is also a
        # cellranger_count output so the DAG waits when both stages run.
        barcodes=all_cellranger_barcodes(config),
    output:
        cells=CC_OUT,
    params:
        script=CC_SCRIPT,
        bam_root=BAM_ROOT,
        samples=",".join(CC_SAMPLES),
        min_umi=CC_MIN_UMI,
        min_genes=CC_MIN_GENES,
    threads: 1
    retries: 0
    resources:
        mem_mb=CC_MEM,
        disk_mb=CC_DISK,
        runtime=CC_RUNTIME,
        lsf_queue=CC_QUEUE,
    shell:
        """
        extra=""
        if [[ -n "{params.min_umi}" ]]; then
            extra="$extra --min-umi {params.min_umi}"
        fi
        if [[ -n "{params.min_genes}" ]]; then
            extra="$extra --min-genes {params.min_genes}"
        fi
        python {params.script} \
            --bam-root {params.bam_root} \
            --samples "{params.samples}" \
            --output {output.cells} \
            $extra
        """
