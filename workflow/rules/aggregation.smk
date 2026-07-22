# Phase 4: tag, combine, and complexity-filter assembled transcripts.

rule edit_cell_transcripts:
    input:
        transcripts="results/cells/{cell}/transcripts.fasta",
    output:
        edited="results/cells/{cell}/transcripts.edit.fasta",
    wildcard_constraints:
        cell="[^/]+",
    shell:
        """
        # Match legacy editing_contigs_array.sh: append _CELL to headers.
        sed "s/>.*/&_{wildcards.cell}/" {input.transcripts} > {output.edited}
        """


rule combine_transcripts:
    input:
        edited=lambda wildcards: all_cell_transcripts_edit(CANDIDATE_CELLS),
    output:
        combined="results/assembly/all_cells.transcripts.edit.fasta",
    shell:
        """
        mkdir -p results/assembly
        cat {input.edited} > {output.combined}
        """


rule filter_low_complexity:
    input:
        combined="results/assembly/all_cells.transcripts.edit.fasta",
    output:
        final=config["outputs"]["combined_transcripts"],
        low_complexity="results/assembly/low_complexity.fq",
    params:
        script="workflow/scripts/run_bbduk_filter.sh",
        entropy=config["params"]["aggregation"]["entropy"],
    threads: config["params"]["aggregation"]["threads"]
    resources:
        mem_mb=config["resources"]["aggregation"]["mem_mb"],
        disk_mb=config["resources"]["aggregation"]["disk_mb"],
        runtime=config["resources"]["aggregation"]["runtime"],
        lsf_queue=config["resources"]["aggregation"]["lsf_queue"],
    shell:
        """
        bash {params.script} \
            {input.combined} \
            {output.final} \
            {output.low_complexity} \
            {params.entropy}
        """
