# Host cell annotation from the PR2 and metaPR2 BLAST tables.

rule summarize_blast_18s:
    input:
        pr2=config["outputs"]["blast_pr2"],
        metapr2=config["outputs"]["blast_metapr2"],
    output:
        summary=config["outputs"]["cell_annotation_blast"],
    params:
        script=vihop_script("workflow/scripts/run_blast_summary.sh"),
        replacement=config["taxonomy"]["replacement_file"],
        identity=config["params"]["host"]["blast_identity"],
        length=config["params"]["host"]["blast_length"],
        level=config["taxonomy"]["level"],
        pre_pident=config["params"]["host"]["blast_prefilter_pident"],
        pre_length=config["params"]["host"]["blast_prefilter_length"],
        force_annotation=config["params"]["host"].get("force_annotation", True),
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
            {params.pre_length} \
            {params.force_annotation}
        """


rule summarize_cell_annotations:
    input:
        blast=config["outputs"]["cell_annotation_blast"],
        cells=config["outputs"]["candidate_cells"],
    output:
        full=config["outputs"]["cell_annotation_full"],
        summarized=config["outputs"]["cell_annotation_final"],
    params:
        script=vihop_script("workflow/scripts/run_summarize_annotations.sh"),
        level=config["taxonomy"]["level"],
        percent=config["params"]["host"]["chimeric_percent"],
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
            {params.percent}
        """


# Local. Runs after the host annotation summary.
rule host_18s_qc:
    input:
        blast=config["outputs"]["cell_annotation_blast"],
        summarized=config["outputs"]["cell_annotation_final"],
        candidates=config["outputs"]["candidate_cells"],
    output:
        report="results/host/qc/host_18s_qc.txt",
        html="results/host/qc/host_18s_qc.html",
        overview="results/host/qc/host_18s_qc_overview.tsv",
        taxa="results/host/qc/host_18s_qc_taxa.tsv",
        pairs="results/host/qc/host_18s_qc_mixed_pairs.tsv",
    params:
        script=vihop_script("workflow/scripts/host_18s_qc.py"),
        outdir="results/host/qc",
        level=config["taxonomy"]["level"],
        percent=config["params"]["host"]["chimeric_percent"],
    shell:
        """
        python {params.script} \
            --blast-summary {input.blast} \
            --summarized {input.summarized} \
            --candidates {input.candidates} \
            --outdir {params.outdir} \
            --level {params.level} \
            --chimeric-percent {params.percent}
        """
