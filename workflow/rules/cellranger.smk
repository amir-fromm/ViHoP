# Cell Ranger count: one job per sample from docs/sample_table_fastq.tsv.
# Multi-run samples use comma-separated FASTQ dirs in column 2 → --fastqs=a,b
#
# LSF (snakemake-executor-plugin-lsf): rusage mem = resources.mem_mb / threads
# (mem_mb is TOTAL job memory). Example huge-ref: threads=8, mem_mb=160000
# → 20 GB/slot × 8 = 160 GB, with --localcores=8 --localmem=160.

from common import (
    all_cellranger_bams,
    cellranger_bam,
    cellranger_barcodes,
    resolve_bam_root,
    resolve_cellranger_samples,
    sample_fastqs_csv,
)

# Absolute script path so jobs work when --directory / --workdir ≠ repo.
CR_SCRIPT = vihop_script("workflow/scripts/run_cellranger_count.sh")

CR_SAMPLES = resolve_cellranger_samples(config)
CR_CFG = config["params"]["cellranger"]
CR_RES = config["resources"]["cellranger"]
CR_TRANSCRIPTOME = config["outputs"]["cellranger_ref"]
# Default bam_root = results/cellranger; Cell Ranger writes <bam_root>/<sample>/
CR_OUTDIR = resolve_bam_root(config)

# Optional advanced knobs (None / null / false → omit from CLI except create_bam).
CR_CREATE_BAM = str(CR_CFG.get("create_bam", True)).lower()
CR_CHEMISTRY = CR_CFG.get("chemistry", "auto")
CR_EXPECT = CR_CFG.get("expect_cells", None)
CR_EXPECT_STR = "" if CR_EXPECT in (None, "null", False, "false") else str(CR_EXPECT)
CR_CHEM_STR = "" if CR_CHEMISTRY in (None, "null", False, "false") else str(CR_CHEMISTRY)

CR_MODULE = CR_CFG.get("module", "") or ""
CR_LOCALMEM = CR_CFG.get("localmem_gb", None)
CR_LOCALMEM_STR = (
    ""
    if CR_LOCALMEM in (None, "null", False, "false")
    else str(int(CR_LOCALMEM))
)
CR_LOCALVMEM = CR_CFG.get("localvmem_gb", None)
CR_LOCALVMEM_STR = (
    ""
    if CR_LOCALVMEM in (None, "null", False, "false")
    else str(int(CR_LOCALVMEM))
)
CR_OVERRIDES = CR_CFG.get("overrides_json", "") or ""
CR_FORCE_CLEAN = str(CR_CFG.get("force_clean", False)).lower()

# Output paths must be strings (not callables) — use {sample} wildcard.
CR_BAM_OUT = f"{CR_OUTDIR}/{{sample}}/outs/possorted_genome_bam.bam"
CR_BARCODES_OUT = (
    f"{CR_OUTDIR}/{{sample}}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
)


rule cellranger_count:
    input:
        ref_done=config["outputs"]["cellranger_ref_done"],
    output:
        # Do not wrap in protected(): breaks Snakemake 9 LSF --target-jobs
        # (AttributeError wildcards is None). keep-incomplete in the LSF
        # profile retains the BAM if the job fails after writing it.
        bam=CR_BAM_OUT,
        barcodes=CR_BARCODES_OUT,
    params:
        script=CR_SCRIPT,
        transcriptome=CR_TRANSCRIPTOME,
        fastqs=lambda wildcards: sample_fastqs_csv(config, wildcards.sample),
        outdir=CR_OUTDIR,
        create_bam=CR_CREATE_BAM,
        chemistry=CR_CHEM_STR,
        expect_cells=CR_EXPECT_STR,
        module=CR_MODULE,
        localmem_gb=CR_LOCALMEM_STR,
        localvmem_gb=CR_LOCALVMEM_STR,
        overrides_json=CR_OVERRIDES,
        force_clean=CR_FORCE_CLEAN,
        mem_mb=CR_RES["mem_mb"],
    threads: CR_CFG["threads"]
    retries: 0
    resources:
        mem_mb=CR_RES["mem_mb"],
        disk_mb=CR_RES["disk_mb"],
        runtime=CR_RES["runtime"],
        lsf_queue=CR_RES["lsf_queue"],
    wildcard_constraints:
        sample="|".join(CR_SAMPLES) if CR_SAMPLES else "no_samples",
    shell:
        """
        bash {params.script} \
            {wildcards.sample} \
            {params.transcriptome} \
            "{params.fastqs}" \
            {params.outdir} \
            {threads} \
            {params.mem_mb} \
            {params.create_bam} \
            "{params.chemistry}" \
            "{params.expect_cells}" \
            "{params.module}" \
            "{params.localmem_gb}" \
            "{params.localvmem_gb}" \
            "{params.overrides_json}" \
            {params.force_clean}
        """


rule cellranger_all:
    input:
        all_cellranger_bams(config),
