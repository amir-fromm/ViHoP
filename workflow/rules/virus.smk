# Virus contig search and host–virus pairing.

rule virus_lastal:
    input:
        transcripts=config["outputs"]["combined_transcripts"],
    output:
        summary=config["outputs"]["viral_contigs"],
    params:
        script="workflow/scripts/run_virus_lastal.sh",
        virus_db=config["references"]["virus_last_db"],
        general_db=config["references"]["general_last_db"],
        workdir="results/virus",
    threads: config["params"]["virus"]["threads"]
    resources:
        mem_mb=config["resources"]["virus"]["mem_mb"],
        disk_mb=config["resources"]["virus"]["disk_mb"],
        runtime=config["resources"]["virus"]["runtime"],
        lsf_queue=config["resources"]["virus"]["lsf_queue"],
    shell:
        """
        bash {params.script} \
            {input.transcripts} \
            {params.virus_db} \
            {params.general_db} \
            {output.summary} \
            {threads} \
            {params.workdir}
        """


rule pair_host_virus:
    input:
        virus=config["outputs"]["viral_contigs"],
        annotation=config["outputs"]["cell_annotation_final"],
    output:
        pairs=config["outputs"]["host_virus_pairs"],
    params:
        script="workflow/scripts/run_virus_pairing.sh",
        taxa=config["references"]["virus_taxonomy"],
        out_dir="results/pairing",
        out_name=config["outputs"]["host_virus_pairs"].rsplit("/", 1)[-1],
        identity=config["params"]["virus"]["identity"],
        length=config["params"]["virus"]["length"],
    resources:
        mem_mb=config["resources"]["virus"]["summary_mem_mb"],
        disk_mb=config["resources"]["virus"]["disk_mb"],
        runtime=120,
        lsf_queue=config["resources"]["virus"]["lsf_queue"],
    shell:
        """
        bash {params.script} \
            {input.virus} \
            {input.annotation} \
            {params.taxa} \
            {params.out_dir} \
            {params.out_name} \
            {params.identity} \
            {params.length}
        """
