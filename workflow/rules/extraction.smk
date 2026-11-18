_SPLIT_SCRIPT = vihop_script("workflow/scripts/split_candidate_cells.py")


rule split_candidate_cells:
    input:
        candidates=config["outputs"]["candidate_cells"],
    output:
        expand("results/candidate_cells_{sample}.tsv", sample=SAMPLES),
    params:
        script=_SPLIT_SCRIPT,
    shell:
        """
        python {params.script} \
            {input.candidates} {output}
        """
