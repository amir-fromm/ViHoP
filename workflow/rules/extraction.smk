rule split_candidate_cells:
    input:
        candidates=config["outputs"]["candidate_cells"],
    output:
        expand("results/candidate_cells_{sample}.tsv", sample=SAMPLES),
    shell:
        """
        python workflow/scripts/split_candidate_cells.py \
            {input.candidates} {output}
        """
