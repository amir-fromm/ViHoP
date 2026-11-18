# Transcriptome reference build: dedup (cdhit|mmseqs|none) → clean → GTF → cellranger mkref.
# Hoist config lookups to module scope — nested config["params"][...] inside a
# params: block confuses the Snakemake DSL parser ("Unexpected keyword config").

REF_NAME = config["references"]["cellranger_ref_name"]
REF_DIR = f"results/ref/{REF_NAME}"
REF_CFG = config["params"]["reference"]
REF_RES = config["resources"]["reference"]
DEDUP_METHOD = REF_CFG["dedup_method"]
DEDUP_IDENTITY = REF_CFG.get("dedup_identity", 0.95)
MMSEQS_CFG = REF_CFG.get("mmseqs", {}) or {}
MMSEQS_COV = MMSEQS_CFG.get("coverage", 0.30)
MMSEQS_COV_MODE = MMSEQS_CFG.get("cov_mode", 1)
MMSEQS_MIN_ALN = MMSEQS_CFG.get("min_aln_len", 200)
MMSEQS_WORKDIR = f"{REF_DIR}/mmseqs_work"
if DEDUP_METHOD == "mmseqs":
    DEDUP_THREADS = int(MMSEQS_CFG.get("threads", 32))
elif DEDUP_METHOD == "cdhit":
    DEDUP_THREADS = int((REF_CFG.get("cdhit") or {}).get("threads", 48))
else:
    DEDUP_THREADS = 1
MKREF_THREADS = REF_CFG["mkref_threads"]
TRANSCRIPTOME_FASTA = config["references"]["transcriptome_fasta"]
MKREF_DONE = config["outputs"]["cellranger_ref_done"]
MKREF_DIR = config["outputs"]["cellranger_ref"]

# none / skip / copy = already-deduplicated FASTA (symlink or copy into place).
if DEDUP_METHOD not in ("cdhit", "mmseqs", "none", "skip", "copy"):
    raise ValueError(
        f"params.reference.dedup_method must be 'cdhit', 'mmseqs', or 'none' "
        f"(skip dedup), got {DEDUP_METHOD!r}"
    )


rule dedup_transcriptome:
    input:
        fasta=TRANSCRIPTOME_FASTA,
    output:
        dedup=f"{REF_DIR}/dedup.fasta",
    params:
        method=DEDUP_METHOD,
        identity=DEDUP_IDENTITY,
        mmseqs_cov=MMSEQS_COV,
        mmseqs_cov_mode=MMSEQS_COV_MODE,
        mmseqs_min_aln=MMSEQS_MIN_ALN,
        mmseqs_workdir=MMSEQS_WORKDIR,
        cdhit_script=vihop_script("workflow/scripts/run_cdhit_dedup.sh"),
        mmseqs_script=vihop_script("workflow/scripts/run_mmseqs_dedup.sh"),
    threads: DEDUP_THREADS
    resources:
        mem_mb=REF_RES["dedup_mem_mb"],
        disk_mb=REF_RES["disk_mb"],
        runtime=REF_RES["runtime"],
        lsf_queue=REF_RES["lsf_queue"],
    shell:
        """
        mkdir -p "$(dirname {output.dedup})"
        if [[ "{params.method}" == "none" || "{params.method}" == "skip" || "{params.method}" == "copy" ]]; then
            echo "[dedup] skip clustering — using input as-is"
            # Prefer hard link / symlink; fall back to copy for cross-filesystem.
            ln -f {input.fasta} {output.dedup} 2>/dev/null \
                || ln -sf "$(cd "$(dirname {input.fasta})" && pwd)/$(basename {input.fasta})" {output.dedup} \
                || cp -f {input.fasta} {output.dedup}
        elif [[ "{params.method}" == "mmseqs" ]]; then
            bash {params.mmseqs_script} \
                {input.fasta} \
                {output.dedup} \
                {params.identity} \
                {threads} \
                {params.mmseqs_cov} \
                {params.mmseqs_cov_mode} \
                {params.mmseqs_min_aln} \
                {params.mmseqs_workdir}
        else
            bash {params.cdhit_script} \
                {input.fasta} \
                {output.dedup} \
                {params.identity} \
                {threads}
        fi
        """


rule clean_transcriptome_fasta:
    input:
        dedup=f"{REF_DIR}/dedup.fasta",
    output:
        final=f"{REF_DIR}/final.fasta",
    params:
        script=vihop_script("workflow/scripts/run_clean_fasta.sh"),
    resources:
        mem_mb=8000,
        disk_mb=REF_RES["disk_mb"],
        runtime=120,
        lsf_queue=REF_RES["lsf_queue"],
    shell:
        """
        bash {params.script} {input.dedup} {output.final}
        """


rule make_pseudo_gtf:
    input:
        fasta=f"{REF_DIR}/final.fasta",
    output:
        gtf=f"{REF_DIR}/genes.gtf",
    params:
        script=vihop_script("workflow/scripts/run_create_gtf.sh"),
    resources:
        mem_mb=8000,
        disk_mb=REF_RES["disk_mb"],
        runtime=240,
        lsf_queue=REF_RES["lsf_queue"],
    shell:
        """
        bash {params.script} {input.fasta} {output.gtf}
        """


rule cellranger_mkref:
    input:
        fasta=f"{REF_DIR}/final.fasta",
        gtf=f"{REF_DIR}/genes.gtf",
    output:
        # Only a done flag — do not declare directory() outputs; Cell Ranger
        # writes mkref_<name>/ first and Snakemake retries + SIGINT then flap.
        done=MKREF_DONE,
    params:
        script=vihop_script("workflow/scripts/run_cellranger_mkref.sh"),
        ref_name=REF_NAME,
        outdir=REF_DIR,
        mem_mb=REF_RES["mkref_mem_mb"],
        # Expected final transcriptome path (for logging / consumers).
        ref_dir=MKREF_DIR,
    threads: MKREF_THREADS
    # Do not auto-retry: a mid-run bkill + retry looks like flapping and
    # leaves orphaned Martian stages.
    retries: 0
    resources:
        mem_mb=REF_RES["mkref_mem_mb"],
        disk_mb=REF_RES["disk_mb"],
        runtime=REF_RES["runtime"],
        lsf_queue=REF_RES["lsf_queue"],
    shell:
        """
        # Resolve absolute paths so mkref (run from outdir) finds inputs.
        fasta="$(cd "$(dirname {input.fasta})" && pwd)/$(basename {input.fasta})"
        gtf="$(cd "$(dirname {input.gtf})" && pwd)/$(basename {input.gtf})"
        bash {params.script} \
            {params.ref_name} \
            "$fasta" \
            "$gtf" \
            {params.outdir} \
            {threads} \
            {params.mem_mb} \
            {output.done}
        """
