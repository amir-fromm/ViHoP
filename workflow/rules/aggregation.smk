# Phase 4: combine edited transcripts (headers are tagged inside assemble_batch).

rule combine_transcripts:
    input:
        batches=expand("results/assembly/batches/{batch}.done", batch=ASSEMBLY_BATCHES),
        candidates=CANDIDATE_CELLS,
    output:
        combined="results/assembly/all_cells.transcripts.edit.fasta",
    params:
        script=vihop_script("workflow/scripts/combine_edited_transcripts.py"),
    shell:
        """
        python {params.script} {input.candidates} {output.combined}
        """


rule filter_low_complexity:
    input:
        combined="results/assembly/all_cells.transcripts.edit.fasta",
    output:
        final=config["outputs"]["combined_transcripts"],
        low_complexity="results/assembly/low_complexity.fq",
    params:
        script=vihop_script("workflow/scripts/run_bbduk_filter.sh"),
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
            {params.entropy} \
            {resources.mem_mb}
        """


# Local. Depends on the filtered FASTA so it runs after BBDuk.
rule assembly_qc:
    input:
        final=config["outputs"]["combined_transcripts"],
        low_complexity="results/assembly/low_complexity.fq",
        candidates=CANDIDATE_CELLS,
    output:
        tsv="results/assembly/assembly_qc.tsv",
        txt="results/assembly/assembly_qc.txt",
    params:
        script=vihop_script("workflow/scripts/assembly_qc.py"),
        cells_dir=config["outputs"]["cells_dir"],
    threads: 4
    shell:
        """
        python {params.script} \
            --fasta {input.final} \
            --candidates {input.candidates} \
            --cells-dir {params.cells_dir} \
            --low-complexity {input.low_complexity} \
            --output-tsv {output.tsv} \
            --output-txt {output.txt} \
            --threads {threads}
        """
