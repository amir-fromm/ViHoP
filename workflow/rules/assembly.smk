# Per-cell RNA-SPAdes assembly (Option A: single reads.fq.gz via --s1).
rule assemble_cell:
    input:
        reads="results/cells/{cell}/reads.fq.gz",
    output:
        transcripts="results/cells/{cell}/transcripts.fasta",
    params:
        script="workflow/scripts/run_rnaspades_cell.sh",
    threads: config["params"]["assembly"]["threads"]
    resources:
        mem_mb=config["resources"]["assembly"]["mem_mb"],
        disk_mb=config["resources"]["assembly"]["disk_mb"],
        runtime=config["resources"]["assembly"]["runtime"],
        lsf_queue=config["resources"]["assembly"]["lsf_queue"],
    wildcard_constraints:
        cell="[^/]+",
    shell:
        """
        bash {params.script} {input.reads} {output.transcripts} {threads}
        """
