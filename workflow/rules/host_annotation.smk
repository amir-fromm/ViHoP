# Host 18S annotation: Barrnap, then BLASTN against PR2 and metaPR2.
# contig_length applies only here. Virus DIAMOND reads the combined
# transcripts, including contigs shorter than this cutoff.

rule filter_transcript_length:
    input:
        fasta=config["outputs"]["combined_transcripts"],
    output:
        filtered="results/host/all_cells.transcripts.filtered.fasta",
    params:
        min_len=config["params"]["host"]["contig_length"],
    shell:
        """
        mkdir -p results/host
        if command -v seqkit >/dev/null 2>&1; then
            seqkit seq -m {params.min_len} {input.fasta} -o {output.filtered}
        else
            # Fallback without seqkit: copy all (length filter skipped with warning).
            echo "WARNING: seqkit not found; copying unfiltered transcripts" >&2
            cp {input.fasta} {output.filtered}
        fi
        """


# One split feeds two parallel searches: Barrnap+BLAST and DIAMOND.
# cmsearch does not split work across nodes, and one DIAMOND blastx does not either.
# params.host.split_fasta (preferred); barrnap_batches kept as an alias.
_HOST_PARAMS = config["params"]["host"]
_BARRNAP_N = int(
    _HOST_PARAMS.get("split_fasta")
    or _HOST_PARAMS.get("barrnap_batches")
    or 40
)
BARRNAP_BATCHES = [f"batch_{i:04d}" for i in range(_BARRNAP_N)]
# cmsearch RSS grows with --cpu. A MetaT batch peaked near 12 GB at 4 threads;
# 64 threads under 16 GB was killed (TERM_MEMLIMIT). Reserve at least
# threads * barrnap_mem_mb_per_thread (default 3000), and never less than
# barrnap_mem_mb.
_BARRNAP_THREADS = max(1, int(config["params"]["host"].get("barrnap_threads", 4)))
_BARRNAP_MEM_FLOOR = int(config["resources"]["host"].get("barrnap_mem_mb", 16000))
_BARRNAP_MEM_PER_THREAD = int(
    config["resources"]["host"].get("barrnap_mem_mb_per_thread", 3000)
)
_BARRNAP_MEM_MB = max(_BARRNAP_MEM_FLOOR, _BARRNAP_THREADS * _BARRNAP_MEM_PER_THREAD)


rule split_barrnap_fasta:
    input:
        # ancient: a later touch of the filtered FASTA must not rewrite every
        # batch and force Barrnap/BLAST to start over.
        filtered=ancient("results/host/all_cells.transcripts.filtered.fasta"),
    output:
        # Not temp(): mid-run restarts must keep the splits.
        fastas=expand("results/host/barrnap_batches/{batch}.fasta", batch=BARRNAP_BATCHES),
    params:
        script=vihop_script("workflow/scripts/split_fasta_batches.py"),
        keep=vihop_script("workflow/scripts/snakemake_keep.sh"),
        keep_dir="results/host/barrnap_batches/_keep",
        outdir="results/host/barrnap_batches",
        n=_BARRNAP_N,
    shell:
        """
        # Snakemake deletes outputs before a re-run. Restore from _keep first.
        bash {params.keep} restore {params.keep_dir} {output.fastas}
        if bash {params.keep} all_present {output.fastas}; then
            echo "[split] keeping existing barrnap batches"
            bash {params.keep} save {params.keep_dir} {output.fastas}
            exit 0
        fi
        python {params.script} {input.filtered} {params.outdir} {params.n}
        bash {params.keep} save {params.keep_dir} {output.fastas}
        """


rule barrnap_batch:
    input:
        fasta="results/host/barrnap_batches/{batch}.fasta",
    output:
        # Not temp(): finished parts must survive a controller restart.
        rrna_fa="results/host/barrnap_parts/{batch}.rrna.fa",
        rrna_gff="results/host/barrnap_parts/{batch}.rrna.gff",
        rrna_18s="results/host/barrnap_parts/{batch}.18s.fasta",
    params:
        script=vihop_script("workflow/scripts/run_barrnap_18s.sh"),
        keep=vihop_script("workflow/scripts/snakemake_keep.sh"),
        keep_dir="results/host/barrnap_parts/_keep",
        kingdom=config["params"]["host"]["barrnap_kingdom"],
    wildcard_constraints:
        batch="batch_[0-9]+",
    threads: _BARRNAP_THREADS
    resources:
        mem_mb=_BARRNAP_MEM_MB,
        disk_mb=config["resources"]["host"]["disk_mb"],
        runtime=config["resources"]["host"]["runtime"],
        lsf_queue=config["resources"]["host"]["lsf_queue"],
    shell:
        """
        # Snakemake deletes these outputs before a re-run; _keep is not a
        # rule output so it survives. Restore, then skip cmsearch if present.
        bash {params.keep} restore {params.keep_dir} \
            {output.rrna_fa} {output.rrna_gff} {output.rrna_18s}
        if bash {params.keep} all_present \
            {output.rrna_fa} {output.rrna_gff} {output.rrna_18s}; then
            echo "[barrnap] keeping existing {output.rrna_18s}"
            bash {params.keep} save {params.keep_dir} \
                {output.rrna_fa} {output.rrna_gff} {output.rrna_18s}
            exit 0
        fi
        bash {params.script} \
            {input.fasta} \
            {output.rrna_fa} \
            {output.rrna_gff} \
            {output.rrna_18s} \
            {threads} \
            {params.kingdom}
        bash {params.keep} save {params.keep_dir} \
            {output.rrna_fa} {output.rrna_gff} {output.rrna_18s}
        """


rule blast_18s_pr2:
    input:
        query="results/host/barrnap_parts/{batch}.18s.fasta",
    output:
        # Not temp(): keep PR2 parts until the merge reads them.
        tsv="results/host/blast_parts/{batch}.PR2.tsv",
    params:
        script=vihop_script("workflow/scripts/run_blastn.sh"),
        keep=vihop_script("workflow/scripts/snakemake_keep.sh"),
        keep_dir="results/host/blast_parts/_keep",
        db=config["references"]["pr2_db"],
        evalue=config["params"]["host"]["blast_evalue"],
    wildcard_constraints:
        batch="batch_[0-9]+",
    threads: config["params"]["host"]["blast_threads"]
    resources:
        mem_mb=config["resources"]["host"]["blast_mem_mb"],
        disk_mb=config["resources"]["host"]["disk_mb"],
        runtime=config["resources"]["host"]["runtime"],
        lsf_queue=config["resources"]["host"]["lsf_queue"],
    shell:
        """
        bash {params.keep} restore {params.keep_dir} {output.tsv}
        if bash {params.keep} all_present {output.tsv}; then
            echo "[blastn] keeping existing {output.tsv}"
            bash {params.keep} save {params.keep_dir} {output.tsv}
            exit 0
        fi
        bash {params.script} \
            {input.query} \
            {params.db} \
            {output.tsv} \
            {params.evalue} \
            {threads}
        bash {params.keep} save {params.keep_dir} {output.tsv}
        """


# Every 18S sequence is searched in metaPR2, including sequences that already
# hit PR2. This uses the same per-batch 18S as PR2, so the two searches run
# together. A tied PR2 hit is dropped later, and only then is metaPR2 needed.


# results/finished/steps.yaml records barrnap_blast after a real merge.
# While that entry is true and the merged tables exist, do not rebuild the
# temporary Barrnap/BLAST batches. Set barrnap_blast: false to run them again.
from common import load_finished_steps

_FINISHED_STEPS = load_finished_steps()
_HOST_BLAST = [
    config["outputs"]["blast_pr2"],
]
_HOST_FILTERED = [
    "results/host/filtered.rRNA.PR2.tsv",
    "results/host/filtered.rRNA.metaPR2.tsv",
]
_SKIP_BARRNAP_BLAST = _FINISHED_STEPS.get("barrnap_blast", False) and (
    all(Path(path).is_file() and Path(path).stat().st_size > 0 for path in _HOST_BLAST)
    or all(Path(path).is_file() and Path(path).stat().st_size > 0 for path in _HOST_FILTERED)
)
# Merged products are the data. When they are already on disk, no rule owns
# them, so a later rerun cannot delete them before the summary.
_MERGED_HOST_FILES = _HOST_BLAST + [
    "results/host/rrna.fa",
    "results/host/rrna.gff",
    config["outputs"]["rrna_18s"],
]
_MERGED_HOST_PRESENT = _FINISHED_STEPS.get("barrnap_blast", False) and all(
    Path(path).is_file() and Path(path).stat().st_size > 0 for path in _MERGED_HOST_FILES
)


def _merge_18s_inputs(wildcards):
    if _SKIP_BARRNAP_BLAST:
        return {"done": ancient("results/finished/steps.yaml")}
    return {
        "pr2": expand("results/host/blast_parts/{batch}.PR2.tsv", batch=BARRNAP_BATCHES),
        "rrna_fa": expand("results/host/barrnap_parts/{batch}.rrna.fa", batch=BARRNAP_BATCHES),
        "rrna_gff": expand("results/host/barrnap_parts/{batch}.rrna.gff", batch=BARRNAP_BATCHES),
        "rrna_18s": expand("results/host/barrnap_parts/{batch}.18s.fasta", batch=BARRNAP_BATCHES),
    }


# Concatenate only after every batch has been searched.
# Omitted once the merged files exist: they are inputs, not a job to rerun.
if not _MERGED_HOST_PRESENT:
  rule merge_18s_blast:
    input:
        unpack(_merge_18s_inputs),
    output:
        # Local rule, so protected() does not hit the LSF --target-jobs bug.
        # Snakemake must not delete these if a later step has not finished.
        # Batch parts are kept on disk (not temp) so a restart can resume.
        pr2=protected(config["outputs"]["blast_pr2"]),
        rrna_fa=protected("results/host/rrna.fa"),
        rrna_gff=protected("results/host/rrna.gff"),
        rrna_18s=protected(config["outputs"]["rrna_18s"]),
    params:
        mark=vihop_script("workflow/scripts/mark_step_finished.py"),
    shell:
        (
            """
            echo "[host] barrnap+blast already finished; restoring merged tables"
            if [[ ! -s {output.pr2} ]]; then
                cp -f results/host/filtered.rRNA.PR2.tsv {output.pr2}
            fi
            if [[ ! -s {output.rrna_fa} ]]; then
                echo "placeholder" > {output.rrna_fa}
            fi
            if [[ ! -s {output.rrna_gff} ]]; then
                echo "placeholder" > {output.rrna_gff}
            fi
            if [[ ! -s {output.rrna_18s} ]]; then
                echo "placeholder" > {output.rrna_18s}
            fi
            """
            if _SKIP_BARRNAP_BLAST
            else """
            mkdir -p results/host
            cat {input.pr2} > {output.pr2}
            cat {input.rrna_fa} > {output.rrna_fa}
            cat {input.rrna_18s} > {output.rrna_18s}
            : > {output.rrna_gff}
            first=1
            for gff in {input.rrna_gff}; do
                if [[ "$first" -eq 1 ]]; then
                    cat "$gff" >> {output.rrna_gff}
                    first=0
                else
                    grep -v '^#' "$gff" >> {output.rrna_gff} || true
                fi
            done
            python {params.mark} barrnap_blast
            """
        )

  # Same 18S batches as PR2. Neither search waits for the other.
  rule blast_18s_metapr2:
    input:
        query="results/host/barrnap_parts/{batch}.18s.fasta",
    output:
        tsv="results/host/blast_parts/{batch}.metaPR2.tsv",
    params:
        script=vihop_script("workflow/scripts/run_blastn.sh"),
        keep=vihop_script("workflow/scripts/snakemake_keep.sh"),
        keep_dir="results/host/blast_parts/_keep",
        db=config["references"]["metapr2_db"],
        evalue=config["params"]["host"]["blast_evalue"],
    wildcard_constraints:
        batch="batch_[0-9]+",
    threads: config["params"]["host"]["blast_threads"]
    resources:
        mem_mb=config["resources"]["host"]["blast_mem_mb"],
        disk_mb=config["resources"]["host"]["disk_mb"],
        runtime=config["resources"]["host"]["runtime"],
        lsf_queue=config["resources"]["host"]["lsf_queue"],
    shell:
        """
        bash {params.keep} restore {params.keep_dir} {output.tsv}
        if bash {params.keep} all_present {output.tsv}; then
            echo "[blastn] keeping existing {output.tsv}"
            bash {params.keep} save {params.keep_dir} {output.tsv}
            exit 0
        fi
        bash {params.script} \
            {input.query} \
            {params.db} \
            {output.tsv} \
            {params.evalue} \
            {threads}
        bash {params.keep} save {params.keep_dir} {output.tsv}
        """

  rule merge_metapr2:
    input:
        parts=expand("results/host/blast_parts/{batch}.metaPR2.tsv", batch=BARRNAP_BATCHES),
    output:
        tsv=protected(config["outputs"]["blast_metapr2"]),
    shell:
        """
        mkdir -p results/host
        cat {input.parts} > {output.tsv}
        """
else:
  # Barrnap+PR2 already finished. Only search metaPR2 if the merged table is
  # still missing. If it exists (often write-protected from merge_metapr2), do
  # not claim it as a rule output — otherwise a later ./vihop downstream run
  # hits ProtectedOutputException before pairing.
  _METAPR2_PATH = Path(config["outputs"]["blast_metapr2"])
  _METAPR2_PRESENT = (
      _METAPR2_PATH.is_file() and _METAPR2_PATH.stat().st_size > 0
  )
  if not _METAPR2_PRESENT:
    rule blast_18s_metapr2:
      input:
          query=config["outputs"]["rrna_18s"],
      output:
          tsv=config["outputs"]["blast_metapr2"],
      params:
          script=vihop_script("workflow/scripts/run_blastn.sh"),
          db=config["references"]["metapr2_db"],
          evalue=config["params"]["host"]["blast_evalue"],
      threads: config["params"]["host"]["blast_threads"]
      resources:
          mem_mb=config["resources"]["host"]["blast_mem_mb"],
          disk_mb=config["resources"]["host"]["disk_mb"],
          runtime=config["resources"]["host"]["runtime"],
          lsf_queue=config["resources"]["host"]["lsf_queue"],
      shell:
          """
          bash {params.script} \
              {input.query} \
              {params.db} \
              {output.tsv} \
              {params.evalue} \
              {threads}
          """

