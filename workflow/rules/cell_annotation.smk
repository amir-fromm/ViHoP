# Host cell annotation: BLAST (+ optional LAST) summarization.

rule summarize_blast_18s:
    input:
        pr2=config["outputs"]["blast_pr2"],
        metapr2=config["outputs"]["blast_metapr2"],
    output:
        summary=config["outputs"]["cell_annotation_blast"],
    params:
        script="workflow/scripts/run_blast_summary.sh",
        replacement=config["taxonomy"]["replacement_file"],
        identity=config["params"]["host"]["blast_identity"],
        length=config["params"]["host"]["blast_length"],
        level=config["taxonomy"]["level"],
        pre_pident=config["params"]["host"]["blast_prefilter_pident"],
        pre_length=config["params"]["host"]["blast_prefilter_length"],
    threads: config["params"]["host"]["summary_threads"]
    resources:
        mem_mb=config["resources"]["host"]["summary_mem_mb"],
        disk_mb=config["resources"]["host"]["disk_mb"],
        runtime=config["resources"]["host"]["runtime"],
        lsf_queue=config["resources"]["host"]["lsf_queue"],
    shell:
        """
        bash {params.script} \
            {input.pr2} \
            {input.metapr2} \
            {output.summary} \
            {params.replacement} \
            {params.identity} \
            {params.length} \
            {params.level} \
            {params.pre_pident} \
            {params.pre_length}
        """


rule summarize_last_18s:
    input:
        last=lambda wildcards: all_cell_lastal_18s(CANDIDATE_CELLS),
    output:
        summary=config["outputs"]["cell_annotation_last"],
    params:
        script="bin/last_sum_18s_multi.py",
        last_dir="results/host/lastal/",
        replacement=config["taxonomy"]["replacement_file"],
        identity=config["params"]["host"]["last_identity"],
        length=config["params"]["host"]["last_length"],
        level=config["taxonomy"]["level"],
        last_count=config["params"]["host"]["last_count"],
    threads: config["params"]["host"]["summary_threads"]
    resources:
        mem_mb=config["resources"]["host"]["summary_mem_mb"],
        disk_mb=config["resources"]["host"]["disk_mb"],
        runtime=config["resources"]["host"]["runtime"],
        lsf_queue=config["resources"]["host"]["lsf_queue"],
    shell:
        """
        python {params.script} \
            --last_path {params.last_dir} \
            --last_summary false \
            --replacement_file {params.replacement} \
            --identity {params.identity} \
            --length {params.length} \
            --count {params.last_count} \
            --level {params.level} \
            --workers {threads}
        mv -f last_summary.{params.level}.csv {output.summary}
        """


def _summarize_annotation_inputs(wildcards):
    inputs = {
        "blast": config["outputs"]["cell_annotation_blast"],
        "cells": config["outputs"]["candidate_cells"],
    }
    if RUN_LASTAL:
        inputs["last"] = config["outputs"]["cell_annotation_last"]
    return inputs


rule summarize_cell_annotations:
    input:
        unpack(_summarize_annotation_inputs),
    output:
        full=f"results/cells/full_annotation.{config['taxonomy']['level']}.csv",
        summarized=config["outputs"]["cell_annotation_final"],
    params:
        script="workflow/scripts/run_summarize_annotations.sh",
        level=config["taxonomy"]["level"],
        percent=config["params"]["host"]["chimeric_percent"],
        last_count=config["params"]["host"]["last_count"],
        last_summary=(
            config["outputs"]["cell_annotation_last"]
            if config["params"]["host"].get("run_lastal", False)
            else ""
        ),
    resources:
        mem_mb=2000,
        disk_mb=4000,
        runtime=60,
        lsf_queue=config["resources"]["host"]["lsf_queue"],
    shell:
        """
        bash {params.script} \
            {input.blast} \
            {input.cells} \
            {output.full} \
            {output.summarized} \
            {params.level} \
            {params.percent} \
            "{params.last_summary}" \
            {params.last_count}
        """
