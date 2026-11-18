# Virus search on its own splits of the combined transcripts.
# The 18S contig_length filter does not apply: short contigs are searched.
# Each batch: DIAMOND blastx against the configured marker proteins, then
# DIAMOND blastx of those candidates against RefSeq. No LAST.

from pathlib import Path as _Path

def _virus_marker_stem(fasta: str) -> str:
    name = _Path(fasta).name
    for suffix in (".faa.gz", ".fasta.gz", ".fa.gz", ".faa", ".fasta", ".fa", ".fna"):
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    return _Path(name).stem


def _virus_lsf_cpus_and_mem(cfg: dict) -> tuple[int, int]:
    """Return (cpus, memory_per_slot_mb) for virus DIAMOND batch jobs.

    snakemake-executor-plugin-lsf:
      -n ← threads
      rusage[mem] ← mem_mb_per_cpu (or mem_mb / threads)
      span[ptile=N] only if resources.mpi is set; otherwise span[hosts=1].
    We set mpi + ptile=cpus so -n equals ptile (required on Wexac so the
    job is not split across hosts incorrectly).
    """
    res = (cfg.get("resources") or {}).get("virus") or {}
    params = (cfg.get("params") or {}).get("virus") or {}
    cpus = int(res.get("cpus") or params.get("threads") or 16)
    if cpus < 1:
        raise ValueError("resources.virus.cpus must be >= 1")
    if "memory_per_slot" in res:
        per_slot = int(res["memory_per_slot"])
    elif "mem_mb" in res:
        # Older configs stored the LSF total; convert to the same per-slot budget.
        per_slot = max(1, int(res["mem_mb"]) // cpus)
    else:
        per_slot = 750
    if per_slot < 1:
        raise ValueError("resources.virus.memory_per_slot must be >= 1")
    return cpus, per_slot


_VIRUS_CPUS, _VIRUS_MEM_PER_SLOT = _virus_lsf_cpus_and_mem(config)
_VIRUS_MARKER_STEM = _virus_marker_stem(config["references"]["virus_db_fasta"])
_VIRUS_MARKER_DB = f"results/virus/{_VIRUS_MARKER_STEM}"
_VIRUS_MARKER_DMND = f"{_VIRUS_MARKER_DB}.dmnd"

# A finished virus merge must not block host Barrnap/BLAST. When the merged
# table is already on disk, omit the virus search rules so Snakemake cannot
# raise ProtectedOutputException while planning the host branch.
# Exception: --virus_search / params.stages.virus_search forces the rules into
# the DAG so a focused rebuild can run (vihop clears results/virus first).
_VIRUS_SUMMARY = config["outputs"]["viral_contigs"]
_CLI_STAGES_LIST = [str(s) for s in (config.get("cli_stages") or [])]
_FORCE_VIRUS_SEARCH = (
    "virus_search" in _CLI_STAGES_LIST
    or bool((config.get("params", {}).get("stages") or {}).get("virus_search"))
)
_MERGED_VIRUS_PRESENT = (
    not _FORCE_VIRUS_SEARCH
    and _Path(_VIRUS_SUMMARY).is_file()
    and _Path(_VIRUS_SUMMARY).stat().st_size > 0
)

# Pairing is part of full downstream, not virus_search-only.
if _CLI_STAGES_LIST:
    _RUN_VIRUS_PAIRING = "downstream" in _CLI_STAGES_LIST
else:
    _RUN_VIRUS_PAIRING = bool(
        (config.get("params", {}).get("stages") or {}).get("downstream", True)
    )


if not _MERGED_VIRUS_PRESENT:
  rule split_virus_fasta:
    input:
        # ancient: a later touch of the combined FASTA must not redo DIAMOND.
        # This is the low-complexity-filtered transcriptome, before the 18S
        # contig_length cutoff.
        fasta=ancient(config["outputs"]["combined_transcripts"]),
    output:
        fastas=expand("results/virus/batches/{batch}.fasta", batch=BARRNAP_BATCHES),
    params:
        script=vihop_script("workflow/scripts/split_fasta_batches.py"),
        keep=vihop_script("workflow/scripts/snakemake_keep.sh"),
        keep_dir="results/virus/batches/_keep",
        outdir="results/virus/batches",
        n=_BARRNAP_N,
    shell:
        """
        bash {params.keep} restore {params.keep_dir} {output.fastas}
        if bash {params.keep} all_present {output.fastas}; then
            echo "[split] keeping existing virus batches"
            bash {params.keep} save {params.keep_dir} {output.fastas}
            exit 0
        fi
        python {params.script} {input.fasta} {params.outdir} {params.n}
        bash {params.keep} save {params.keep_dir} {output.fastas}
        """


  rule make_virus_marker_dmnd:
    input:
        fasta=config["references"]["virus_db_fasta"],
    output:
        dmnd=_VIRUS_MARKER_DMND,
    params:
        diamond=config["references"].get("diamond_bin", "diamond"),
        db=_VIRUS_MARKER_DB,
    threads: 4
    resources:
        # -n 4 with span[ptile=4] (plugin needs mpi=True to emit ptile).
        mem_mb_per_cpu=1000,
        mpi=True,
        ptile=4,
        disk_mb=4000,
        runtime=60,
        lsf_queue=config["resources"]["virus"]["lsf_queue"],
    shell:
        """
        mkdir -p results/virus
        {params.diamond} makedb --in "{input.fasta}" --db "{params.db}" --threads {threads}
        """


  rule virus_diamond_batch:
    input:
        transcripts="results/virus/batches/{batch}.fasta",
        marker_dmnd=_VIRUS_MARKER_DMND,
    output:
        # Not temp(): finished parts must survive a controller restart / failed
        # merge. Removed only after viral_contigs is written (merge rule).
        summary="results/virus/parts/{batch}.tsv",
    params:
        script=vihop_script("workflow/scripts/run_virus_diamond.sh"),
        diamond_db=config["references"]["general_diamond_db"],
        catalog=config["references"]["refseq_catalog"],
        nodes=config["references"]["refseq_nodes"],
        merged=config["references"]["refseq_merged"],
        diamond_bin=config["references"].get("diamond_bin", "diamond"),
        workdir="results/virus/work/{batch}",
    wildcard_constraints:
        batch="batch_[0-9]+",
    # Transient LSF EXIT (0s, empty log) is common under high fan-out; retry
    # these batches without touching the global profile retries: 0 (mkref).
    retries: 2
    threads: _VIRUS_CPUS
    resources:
        # LSF: -n cpus, rusage[mem=memory_per_slot], span[ptile=cpus].
        # mpi=True is required by snakemake-executor-plugin-lsf to honor ptile
        # (otherwise it forces span[hosts=1] and never sets ptile=n).
        mem_mb_per_cpu=_VIRUS_MEM_PER_SLOT,
        mpi=True,
        ptile=_VIRUS_CPUS,
        disk_mb=config["resources"]["virus"]["disk_mb"],
        runtime=config["resources"]["virus"]["runtime"],
        lsf_queue=config["resources"]["virus"]["lsf_queue"],
    shell:
        """
        rm -rf {params.workdir}
        bash {params.script} \
            {input.transcripts} \
            {input.marker_dmnd} \
            {params.diamond_db} \
            {output.summary} \
            {threads} \
            {params.workdir} \
            {params.catalog} \
            {params.nodes} \
            {params.merged} \
            {params.diamond_bin}
        rm -rf {params.workdir}
        """


  # Merge only after every split has finished DIAMOND.
  rule merge_virus_batches:
    input:
        parts=expand("results/virus/parts/{batch}.tsv", batch=BARRNAP_BATCHES),
    output:
        summary=protected(config["outputs"]["viral_contigs"]),
    params:
        verify=vihop_script("workflow/scripts/verify_virus_batches.py"),
        mark=vihop_script("workflow/scripts/mark_step_finished.py"),
        n=_BARRNAP_N,
        parts_dir="results/virus/parts",
    shell:
        """
        python {params.verify} {params.n}
        mkdir -p "$(dirname {output.summary})"
        cat {input.parts} > {output.summary}
        python {params.mark} virus_diamond
        # Drop per-batch TSVs only after the final merge exists on disk.
        if [[ -s "{output.summary}" ]]; then
          rm -f {params.parts_dir}/batch_*.tsv
        fi
        """


if _RUN_VIRUS_PAIRING:
  rule pair_host_virus:
    input:
        virus=config["outputs"]["viral_contigs"],
        annotation=config["outputs"]["cell_annotation_final"],
    output:
        pairs=config["outputs"]["host_virus_pairs"],
    params:
        script=vihop_script("workflow/scripts/run_virus_pairing.sh"),
        taxa=config["references"]["virus_taxonomy"],
        out_dir="results/pairing",
        out_name=config["outputs"]["host_virus_pairs"].rsplit("/", 1)[-1],
        identity=config["params"]["virus"]["identity"],
        length=config["params"]["virus"]["length"],
    threads: 1
    resources:
        mem_mb_per_cpu=config["resources"]["virus"]["summary_mem_mb"],
        mpi=True,
        ptile=1,
        disk_mb=config["resources"]["virus"]["disk_mb"],
        runtime=120,
        lsf_queue=config["resources"]["virus"]["lsf_queue"],
    shell:
        """
        bash {params.script} \
            {input.virus} \
            {input.annotation} \
            {params.taxa} \
            {params.out_dir} \
            {params.out_name} \
            {params.identity} \
            {params.length}
        """
