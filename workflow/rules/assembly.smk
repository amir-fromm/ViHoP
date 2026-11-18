# One LSF job assembles many cells. The login node builds this small DAG once;
# the compute job does not rebuild a per-cell graph.
rule assemble_batch:
    input:
        manifest="results/assembly/batches/{batch}.tsv",
        extracts=expand("results/extraction/{sample}.done", sample=SAMPLES),
    output:
        done="results/assembly/batches/{batch}.done",
    params:
        script=vihop_script("workflow/scripts/run_rnaspades_batch.py"),
        cell_script=vihop_script("workflow/scripts/run_rnaspades_cell.sh"),
    wildcard_constraints:
        batch="batch_[0-9]+",
    threads: config["params"]["assembly"]["threads"]
    resources:
        mem_mb=config["resources"]["assembly"]["mem_mb"],
        disk_mb=config["resources"]["assembly"]["disk_mb"],
        runtime=config["resources"]["assembly"]["runtime"],
        lsf_queue=config["resources"]["assembly"]["lsf_queue"],
    shell:
        """
        python {params.script} {input.manifest} {output.done} {threads} {params.cell_script}
        """
