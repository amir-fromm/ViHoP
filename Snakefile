import sys
from pathlib import Path

# Prefer a project-local config (workdir / cwd), then the repo example.
# With --directory, relative paths are resolved against the workdir — so a cloned
# GitHub project can keep its own config.yaml without editing the Snakefile.
_REPO_CFG = Path(workflow.basedir) / "config" / "config.yaml"
_LOCAL_CANDIDATES = (Path("config.yaml"), Path("config/config.yaml"))
_CFG = next((p for p in _LOCAL_CANDIDATES if p.is_file()), None)
if _CFG is None:
    if not _REPO_CFG.is_file():
        raise FileNotFoundError(
            "No config found. Put config.yaml or config/config.yaml in the "
            f"working directory, or keep the example at {_REPO_CFG}"
        )
    _CFG = _REPO_CFG
configfile: str(_CFG)

# Summary filenames follow taxonomy.level (Class, Order, …). A literal
# "Class" in the config is not used.
# ./vihop --annotation --bp N sets annotation_named_outputs and writes
# blast_summary.<level>.<bp>.csv so a second cutoff does not replace the first.
_TAXON_LEVEL = str(config["taxonomy"]["level"])
_HOST_BP = int(config.get("params", {}).get("host", {}).get("blast_length", 300))
if config.get("annotation_named_outputs"):
    _ANN_TAG = f"{_TAXON_LEVEL}.{_HOST_BP}"
else:
    _ANN_TAG = _TAXON_LEVEL
config["outputs"]["cell_annotation_blast"] = f"results/host/blast_summary.{_ANN_TAG}.csv"
config["outputs"]["cell_annotation_full"] = f"results/host/full_annotation.{_ANN_TAG}.csv"
config["outputs"]["cell_annotation_last"] = f"results/host/last_summary.{_TAXON_LEVEL}.csv"
config["outputs"]["cell_annotation_final"] = f"results/host/summarized_annotation.{_ANN_TAG}.csv"
config["outputs"]["host_virus_pairs"] = f"results/pairing/summary_blast_viruses.{_TAXON_LEVEL}.tsv"

# Install root = directory that contains ./vihop and workflow/scripts/.
# Set by the vihop wrapper; fall back to Snakefile location if snakemake is
# invoked directly.
VIHOP_ROOT = Path(config["vihop_root"]) if config.get("vihop_root") else Path(workflow.basedir)
VIHOP_SCRIPTS = VIHOP_ROOT / "workflow" / "scripts"
VIHOP_BIN = VIHOP_ROOT / "bin"


def vihop_script(*parts: str) -> str:
    """Absolute path under VIHOP_ROOT (e.g. vihop_script('workflow/scripts/foo.py'))."""
    return str(VIHOP_ROOT.joinpath(*parts))


sys.path.insert(0, str(VIHOP_SCRIPTS))

from common import apply_lsf_queue_names

apply_lsf_queue_names(config)

STAGES = dict(config.get("params", {}).get("stages", {}) or {})
# ./vihop --choose_cells (etc.): cli_stages list forces only those stages on.
_CLI_STAGES = config.get("cli_stages")
if _CLI_STAGES:
    _allowed = {
        "reference",
        "cellranger",
        "choose_cells",
        "downstream",
        "virus_search",
        "annotation",
    }
    _selected = {str(s) for s in _CLI_STAGES}
    unknown = _selected - _allowed
    if unknown:
        raise ValueError(f"Unknown cli_stages: {sorted(unknown)}; allowed={sorted(_allowed)}")
    for _name in _allowed:
        STAGES[_name] = _name in _selected

RUN_REFERENCE = bool(STAGES.get("reference", False))
RUN_CELLRANGER = bool(STAGES.get("cellranger", False))
RUN_CHOOSE_CELLS = bool(STAGES.get("choose_cells", False))
RUN_DOWNSTREAM = bool(STAGES.get("downstream", True))
RUN_VIRUS_SEARCH = bool(STAGES.get("virus_search", False))
RUN_ANNOTATION = bool(STAGES.get("annotation", False))
# Virus-only: --virus_search without --downstream. With both, full downstream.
RUN_VIRUS_SEARCH_ONLY = bool(RUN_VIRUS_SEARCH and not RUN_DOWNSTREAM)
# Host rerun: --annotation --bp N. Filter, Barrnap, both BLASTs, then summary.
# Does not expand extraction, assembly, virus, or pairing.
RUN_ANNOTATION_ONLY = bool(
    RUN_ANNOTATION
    and not RUN_DOWNSTREAM
    and not RUN_VIRUS_SEARCH
    and not RUN_REFERENCE
    and not RUN_CELLRANGER
    and not RUN_CHOOSE_CELLS
)
# Host annotation is BLAST-only (PR2 + metaPR2). LAST is not used.
RUN_LASTAL = False

# reference.smk only when building a transcriptome; Cell Ranger can use a
# prebuilt ref via outputs.cellranger_ref + cellranger_ref_done (e.g. reference.json).
if RUN_REFERENCE:
    include: "workflow/rules/reference.smk"

if RUN_CELLRANGER:
    include: "workflow/rules/cellranger.smk"

if RUN_CHOOSE_CELLS:
    include: "workflow/rules/choose_cells.smk"

# Downstream expands one extract/assemble rule per row of candidate_cells.tsv at
# parse time. If the file was deleted, rebuild it via choose_cells first, then
# re-run with downstream enabled (file present).
CANDIDATE_CELLS = config["outputs"]["candidate_cells"]
_HAVE_CANDIDATES = Path(CANDIDATE_CELLS).is_file()
RUN_DOWNSTREAM_EXPAND = bool(RUN_DOWNSTREAM and _HAVE_CANDIDATES)

if RUN_DOWNSTREAM and not _HAVE_CANDIDATES:
    if not RUN_CHOOSE_CELLS:
        raise FileNotFoundError(
            f"{CANDIDATE_CELLS} is missing. Enable params.stages.choose_cells "
            "to rebuild it from Cell Ranger filtered matrices, or restore the file."
        )
    print(
        f"NOTE: {CANDIDATE_CELLS} missing — skipping downstream expand this parse "
        "(vihop two-pass runs choose_cells first, then downstream).",
        file=sys.stderr,
    )

if RUN_DOWNSTREAM_EXPAND:
    from common import (
        get_samples,
        pending_lists_ready,
        write_assembly_batches,
        write_pending_cell_lists,
    )
    from generate_extraction_rules import write_extraction_rules

    SAMPLES = get_samples(CANDIDATE_CELLS)
    _BATCH_SIZE = int(config["params"]["assembly"].get("cells_per_batch", 200))
    # ./vihop writes these from results/cells before invoking snakemake and
    # deletes them after it exits. Re-parsing on a cluster job must reuse them.
    if not pending_lists_ready():
        _need_extract, _need_assemble = write_pending_cell_lists(
            CANDIDATE_CELLS,
            config["outputs"]["cells_dir"],
        )
        print(
            f"[pending] wrote tmp_extract ({len(_need_extract)}) "
            f"and tmp_assemble ({len(_need_assemble)})",
            file=sys.stderr,
        )
    ASSEMBLY_BATCHES = write_assembly_batches("results/tmp_assemble.tsv", _BATCH_SIZE)

    EXTRACTION_RULES = VIHOP_ROOT / "workflow/rules/extraction_generated.smk"
    write_extraction_rules(
        EXTRACTION_RULES,
        SAMPLES,
        config,
        "results/tmp_extract.tsv",
    )

    include: "workflow/rules/extraction.smk"
    include: "workflow/rules/extraction_generated.smk"
    include: "workflow/rules/assembly.smk"
    include: "workflow/rules/aggregation.smk"
    include: "workflow/rules/host_annotation.smk"
    include: "workflow/rules/cell_annotation.smk"
    include: "workflow/rules/virus.smk"
elif RUN_VIRUS_SEARCH_ONLY:
    # Virus DIAMOND only. Splits the combined transcripts, not the 18S
    # length-filtered FASTA, and does not plan host Barrnap/BLAST.
    _host_p = config["params"]["host"]
    _BARRNAP_N = int(
        _host_p.get("split_fasta") or _host_p.get("barrnap_batches") or 40
    )
    BARRNAP_BATCHES = [f"batch_{i:04d}" for i in range(_BARRNAP_N)]
    _combined = Path(config["outputs"]["combined_transcripts"])
    if not _combined.is_file() or _combined.stat().st_size == 0:
        raise FileNotFoundError(
            "virus_search requires the combined transcripts at "
            f"{_combined} (no contig_length cutoff). "
            "Finish assembly first. virus_search does not rebuild host Barrnap/BLAST."
        )
    include: "workflow/rules/virus.smk"
elif RUN_ANNOTATION_ONLY:
    if not Path(CANDIDATE_CELLS).is_file():
        raise FileNotFoundError(
            f"{CANDIDATE_CELLS} is missing. Host annotation needs the cell list "
            "written by choose_cells."
        )
    include: "workflow/rules/host_annotation.smk"
    include: "workflow/rules/cell_annotation.smk"

# dedup_transcriptome is local when dedup_method is none/skip/copy (just ln/cp).
_DEDUP_LOCAL = str(
    config.get("params", {}).get("reference", {}).get("dedup_method", "")
).lower() in ("none", "skip", "copy")

if RUN_DOWNSTREAM_EXPAND and RUN_CELLRANGER and RUN_CHOOSE_CELLS:
    localrules:
        all,
        cellranger_all,
        choose_cells,
        split_candidate_cells,
        combine_transcripts,
        assembly_qc,
        filter_transcript_length,
        summarize_cell_annotations,
        host_18s_qc,
        pair_host_virus,
elif RUN_DOWNSTREAM_EXPAND and RUN_CHOOSE_CELLS:
    localrules:
        all,
        choose_cells,
        split_candidate_cells,
        combine_transcripts,
        assembly_qc,
        filter_transcript_length,
        summarize_cell_annotations,
        host_18s_qc,
        pair_host_virus,
elif RUN_DOWNSTREAM_EXPAND and RUN_CELLRANGER:
    localrules:
        all,
        cellranger_all,
        split_candidate_cells,
        combine_transcripts,
        assembly_qc,
        filter_transcript_length,
        summarize_cell_annotations,
        host_18s_qc,
        pair_host_virus,
elif RUN_DOWNSTREAM_EXPAND:
    localrules:
        all,
        split_candidate_cells,
        combine_transcripts,
        assembly_qc,
        filter_transcript_length,
        summarize_cell_annotations,
        host_18s_qc,
        pair_host_virus,
elif RUN_CELLRANGER and RUN_CHOOSE_CELLS:
    localrules:
        all,
        cellranger_all,
        choose_cells,
elif RUN_CHOOSE_CELLS:
    localrules:
        all,
        choose_cells,
elif RUN_CELLRANGER:
    localrules:
        all,
        cellranger_all,
elif RUN_VIRUS_SEARCH_ONLY:
    localrules:
        all,
        split_virus_fasta,
        make_virus_marker_dmnd,
        merge_virus_batches,
elif RUN_ANNOTATION_ONLY:
    localrules:
        all,
        filter_transcript_length,
        split_barrnap_fasta,
        merge_18s_blast,
        merge_metapr2,
        summarize_cell_annotations,
else:
    localrules:
        all,

if RUN_REFERENCE and _DEDUP_LOCAL:
    localrules: dedup_transcriptome

if RUN_DOWNSTREAM_EXPAND and not _MERGED_HOST_PRESENT and not _MERGED_VIRUS_PRESENT:
    localrules:
        split_barrnap_fasta,
        merge_18s_blast,
        merge_metapr2,
        split_virus_fasta,
        make_virus_marker_dmnd,
        merge_virus_batches,
elif RUN_DOWNSTREAM_EXPAND and not _MERGED_HOST_PRESENT:
    localrules:
        split_barrnap_fasta,
        merge_18s_blast,
        merge_metapr2,
elif RUN_DOWNSTREAM_EXPAND and not _MERGED_VIRUS_PRESENT:
    localrules:
        split_barrnap_fasta,
        split_virus_fasta,
        make_virus_marker_dmnd,
        merge_virus_batches,
elif RUN_DOWNSTREAM_EXPAND:
    localrules:
        split_barrnap_fasta,


def final_targets():
    targets = []
    if RUN_REFERENCE and not RUN_CELLRANGER and not RUN_CHOOSE_CELLS:
        # Standalone reference build.
        targets.append(config["outputs"]["cellranger_ref_done"])
    if RUN_CELLRANGER:
        from common import all_cellranger_bams

        # Cell Ranger implies a usable transcriptome (mkref.done as dependency).
        targets.extend(all_cellranger_bams(config))
    if RUN_CHOOSE_CELLS:
        targets.append(config["outputs"]["candidate_cells"])
    if RUN_VIRUS_SEARCH_ONLY:
        targets.append(config["outputs"]["viral_contigs"])
    elif RUN_ANNOTATION_ONLY:
        targets.extend(
            [
                config["outputs"]["cell_annotation_blast"],
                config["outputs"]["cell_annotation_full"],
                config["outputs"]["cell_annotation_final"],
            ]
        )
    elif RUN_DOWNSTREAM_EXPAND:
        targets.extend(
            [
                config["outputs"]["combined_transcripts"],
                "results/assembly/assembly_qc.tsv",
                "results/assembly/assembly_qc.txt",
                config["outputs"]["blast_pr2"],
                config["outputs"]["blast_metapr2"],
                config["outputs"]["cell_annotation_blast"],
                config["outputs"]["cell_annotation_final"],
                "results/host/qc/host_18s_qc.txt",
                "results/host/qc/host_18s_qc.html",
                config["outputs"]["viral_contigs"],
                config["outputs"]["host_virus_pairs"],
            ]
        )
    if not targets:
        raise ValueError(
            "No pipeline targets: set params.stages.reference, "
            "params.stages.cellranger, params.stages.choose_cells, "
            "params.stages.downstream, params.stages.virus_search, "
            "and/or ./vihop --annotation --bp N"
        )
    return targets


rule all:
    input:
        final_targets(),
