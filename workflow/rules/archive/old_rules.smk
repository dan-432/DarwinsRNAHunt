#  Miscelanious archived rules that are no longer used in the main pipeline, but may be useful for future reference.

rule download_interpro_taxonomy:
    output:
        taxonomy = "data/interpro/{target_domain}_taxonomy.json"
    params:
        interpro_id = "{target_domain}",
        interpro_base_url = lambda wildcards: config["databases"]["interpro"]["tax_base_url"],
        script = "workflow/src/download_interpro_taxonomy.py"
    log:
        "logs/{target_domain}_taxonomy/download_interpro_taxonomy_{target_domain}.log"
    conda:
        "envs/domain_taxonomy_env.yaml"
    shell:
        """
        python {params.script} \
            --interpro-id {params.interpro_id} \
            --interpro-tax-base-url {params.interpro_base_url} \
            --output {output.taxonomy} \
            2>&1 | tee {log}
        """