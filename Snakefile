configfile: "config/config.yaml"

import sys
from pathlib import Path

sys.path.insert(0, str(Path(workflow.basedir) / "workflow" / "scripts"))
from common import (
    all_cell_extraction_outputs,
    all_cell_lastal_18s,
    all_cell_transcripts,
    all_cell_transcripts_edit,
    get_samples,
)
from generate_extraction_rules import write_extraction_rules

CANDIDATE_CELLS = config["outputs"]["candidate_cells"]
SAMPLES = get_samples(CANDIDATE_CELLS)
RUN_LASTAL = bool(config["params"]["host"].get("run_lastal", False))

EXTRACTION_RULES = Path(workflow.basedir) / "workflow/rules/extraction_generated.smk"
write_extraction_rules(EXTRACTION_RULES, CANDIDATE_CELLS, SAMPLES, config)

include: "workflow/rules/extraction.smk"
include: "workflow/rules/extraction_generated.smk"
include: "workflow/rules/assembly.smk"
include: "workflow/rules/aggregation.smk"
include: "workflow/rules/host_annotation.smk"
include: "workflow/rules/cell_annotation.smk"
include: "workflow/rules/virus.smk"

localrules:
    all,
    split_candidate_cells,
    edit_cell_transcripts,
    combine_transcripts,
    filter_transcript_length,
    summarize_cell_annotations,
    pair_host_virus,


def final_targets():
    targets = [
        *all_cell_extraction_outputs(CANDIDATE_CELLS),
        *all_cell_transcripts(CANDIDATE_CELLS),
        config["outputs"]["combined_transcripts"],
        config["outputs"]["blast_pr2"],
        config["outputs"]["blast_metapr2"],
        config["outputs"]["cell_annotation_blast"],
        config["outputs"]["cell_annotation_final"],
        config["outputs"]["viral_contigs"],
        config["outputs"]["host_virus_pairs"],
    ]
    if RUN_LASTAL:
        targets.extend(all_cell_lastal_18s(CANDIDATE_CELLS))
        targets.append(config["outputs"]["cell_annotation_last"])
    return targets


rule all:
    input:
        final_targets(),
