# Host 18S annotation searches (BLAST on contigs + LAST on reads).

rule filter_transcript_length:
    input:
        fasta=config["outputs"]["combined_transcripts"],
    output:
        filtered="results/host/all_cells.transcripts.filtered.fasta",
    params:
        min_len=config["params"]["host"]["contig_length"],
    shell:
        """
        mkdir -p results/host
        if command -v seqkit >/dev/null 2>&1; then
            seqkit seq -m {params.min_len} {input.fasta} -o {output.filtered}
        else
            # Fallback without seqkit: copy all (length filter skipped with warning).
            echo "WARNING: seqkit not found; copying unfiltered transcripts" >&2
            cp {input.fasta} {output.filtered}
        fi
        """


rule barrnap_18s:
    input:
        filtered="results/host/all_cells.transcripts.filtered.fasta",
    output:
        rrna_fa="results/host/rrna.fa",
        rrna_gff="results/host/rrna.gff",
        rrna_18s=config["outputs"]["rrna_18s"],
    params:
        script="workflow/scripts/run_barrnap_18s.sh",
        kingdom=config["params"]["host"]["barrnap_kingdom"],
    threads: config["params"]["host"]["barrnap_threads"]
    resources:
        mem_mb=config["resources"]["host"]["barrnap_mem_mb"],
        disk_mb=config["resources"]["host"]["disk_mb"],
        runtime=config["resources"]["host"]["runtime"],
        lsf_queue=config["resources"]["host"]["lsf_queue"],
    shell:
        """
        bash {params.script} \
            {input.filtered} \
            {output.rrna_fa} \
            {output.rrna_gff} \
            {output.rrna_18s} \
            {threads} \
            {params.kingdom}
        """


rule blast_18s_pr2:
    input:
        query=config["outputs"]["rrna_18s"],
    output:
        tsv=config["outputs"]["blast_pr2"],
    params:
        script="workflow/scripts/run_blastn.sh",
        db=config["references"]["pr2_db"],
        evalue=config["params"]["host"]["blast_evalue"],
    threads: config["params"]["host"]["blast_threads"]
    resources:
        mem_mb=config["resources"]["host"]["blast_mem_mb"],
        disk_mb=config["resources"]["host"]["disk_mb"],
        runtime=config["resources"]["host"]["runtime"],
        lsf_queue=config["resources"]["host"]["lsf_queue"],
    shell:
        """
        bash {params.script} \
            {input.query} \
            {params.db} \
            {output.tsv} \
            {params.evalue} \
            {threads}
        """


rule blast_18s_metapr2:
    input:
        query=config["outputs"]["rrna_18s"],
    output:
        tsv=config["outputs"]["blast_metapr2"],
    params:
        script="workflow/scripts/run_blastn.sh",
        db=config["references"]["metapr2_db"],
        evalue=config["params"]["host"]["blast_evalue"],
    threads: config["params"]["host"]["blast_threads"]
    resources:
        mem_mb=config["resources"]["host"]["blast_mem_mb"],
        disk_mb=config["resources"]["host"]["disk_mb"],
        runtime=config["resources"]["host"]["runtime"],
        lsf_queue=config["resources"]["host"]["lsf_queue"],
    shell:
        """
        bash {params.script} \
            {input.query} \
            {params.db} \
            {output.tsv} \
            {params.evalue} \
            {threads}
        """


rule lastal_18s:
    input:
        reads="results/cells/{cell}/reads.fq.gz",
    output:
        last="results/host/lastal/{cell}.last",
    params:
        script="workflow/scripts/run_lastal_18s.sh",
        db=config["references"]["pr2_last_db"],
        identity=config["params"]["host"]["last_identity"],
        length=config["params"]["host"]["last_length"],
        mem_mb=config["resources"]["host"]["last_mem_mb"],
    threads: config["params"]["host"]["last_threads"]
    resources:
        mem_mb=config["resources"]["host"]["last_mem_mb"],
        disk_mb=config["resources"]["host"]["disk_mb"],
        runtime=config["resources"]["host"]["runtime"],
        lsf_queue=config["resources"]["host"]["lsf_queue"],
    wildcard_constraints:
        cell="[^/]+",
    shell:
        """
        bash {params.script} \
            {input.reads} \
            {output.last} \
            {params.db} \
            {params.identity} \
            {params.length} \
            {threads} \
            {params.mem_mb}
        """
