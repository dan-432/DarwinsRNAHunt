"""
Sequence set curation rules for filtering and selecting genomes based on taxonomic and quality criteria.
Extracts flanking sequences around identified protein domains, combines them, and applies filtering based on length and GC content.
"""

def get_accessions(wildcards):
    """Get the list of genome accessions for a given target_domain from the download_genomes checkpoint."""
    from json import load
    # Pass the target_domain wildcard to the checkpoint
    checkpoint_output = checkpoints.download_genomes.get(target_domain=wildcards.target_domain)
    with open(checkpoint_output.output.metadata) as f:
        metadata = load(f)
    return metadata["accessions"]

rule download_gtdb_data:
    output:
        tree = "data/gtdb/bac120.tree",
        metadata = "data/gtdb/bac120_metadata.tsv.gz"
    params:
        tree_url = config["databases"]["gtdb"]["gtdb_bacteria_tree_url"],
        metadata_url = config["databases"]["gtdb"]["gtdb_bacteria_metadata_url"],
        script = "workflow/src/download_gtdb.py"
    log:
        "logs/gtdb/download_gtdb_data.log"
    resources:
        mem_mb=config["genomes"]["download_gtdb"]["mem_mb"],
        disk_mb=config["genomes"]["download_gtdb"]["disk_mb"]
    conda:
        "../envs/domain_taxonomy_env.yaml"
    shell:
        """
        python {params.script} \
            --tree-url {params.tree_url} \
            --metadata-url {params.metadata_url} \
            --tree-output {output.tree} \
            --metadata-output {output.metadata} \
            2>&1 | tee {log}
        """

rule download_rfam:
    output:
        cm = "data/rfam/Rfam.cm",
        genome_index = "data/rfam/Rfam.full_region"
    params:
        rfam_cm_url = config["databases"]["rfam"]["rfam_cm_url"],
        rfam_genome_index_url = config["databases"]["rfam"]["rfam_genome_index_url"],
        script = "workflow/src/download_rfam.py"
    log:
        "logs/rfam/download_rfam.log"
    conda:
        "../envs/domain_taxonomy_env.yaml"
    shell:
        """
        python {params.script} \
            --rfam-cm-url {params.rfam_cm_url} \
            --rfam-genome-index-url {params.rfam_genome_index_url} \
            --output-cm {output.cm} \
            --output-genome-index {output.genome_index} \
            2>&1 | tee {log}
        """

rule extract_rfam_cm:
    input:
        cm = "data/rfam/Rfam.cm"
    output:
        cm = "data/rfam/covariance_models/{family}.cm"
    params:
        key = "{family}"
    conda:
        "../envs/rna_motif_env.yaml"
    log:
        "logs/rfam/extract_{family}_cm.log"
    shell:
        "cmfetch {input.cm} {params.key} > {output.cm} 2> {log}"


rule select_genomes:
    input:
        tree = "data/gtdb/bac120.tree",
        metadata = "data/gtdb/bac120_metadata.tsv.gz"
        #to_keep = "data/interpro/{target_domain}_taxonomy.json"
    output:
        tree = "results/01_species_selection/{target_domain}_taxonomy/filtered_tree.tree",
        metadata = "results/01_species_selection/{target_domain}_taxonomy/filtered_metadata.tsv.gz",
        genomes = "results/01_species_selection/{target_domain}_taxonomy/select_genome_accessions.json"
        # image = "results/taxonomy/filtered_tree.png"
    params:
        script = "workflow/src/extract_subtree.py",
        # Looked up per target_domain, since each domain now carries its own taxon.
        target_order = lambda wildcards: TARGETS[wildcards.target_domain]["taxonomic_node"],
        # Shared across all targets.
        completeness = config["taxonomy"]["completeness_threshold"],
        contamination = config["taxonomy"]["contamination_threshold"]
    log:
        "logs/{target_domain}_taxonomy/select_genomes_{target_domain}.log"
    conda:
        "../envs/domain_taxonomy_env.yaml"
    shell:
        """
        python {params.script} \
            --input-tree {input.tree} \
            --metadata {input.metadata} \
            --target-taxanomic-order "{params.target_order}" \
            --completeness {params.completeness} \
            --contamination {params.contamination} \
            --output-tree {output.tree} \
            --output-metadata {output.metadata} \
            --output-genomes {output.genomes} \
            --from-interpro \
            2>&1 | tee {log}
        """

checkpoint download_genomes:
    input:
        genome_accessions = "results/01_species_selection/{target_domain}_taxonomy/select_genome_accessions.json"
    output:
        metadata = "results/01_species_selection/{target_domain}_taxonomy/genomes_available_metadata.json",
    params:
        script = "workflow/src/download_ncbi_genomes.sh",
        outdir = "data/genomes",
        formats = "gff3,genome,protein",
        batch_size = config["databases"]["ncbi"]["batch_size"],
        genome_log = "data/genomes/genome_log.json",
        temp_dir = ".temp"
    log:
        "logs/genomes/download_genomes_{target_domain}.log"
    conda:
        "../envs/genome_download_env.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        mkdir -p {params.temp_dir}
        bash {params.script} \
            {input.genome_accessions} \
            {params.outdir} \
            {params.genome_log} \
            {output.metadata} \
            {params.formats} \
            {params.batch_size} \
            {params.temp_dir} \
            2>&1 | tee {log}
        """

rule locate_domains:
    input:
        hmm = "data/profiles/{target_domain}_domain_profile.hmm",
        protein = "data/genomes/{genome}/protein.faa",
        flag = lambda wildcards: checkpoints.download_genomes.get(target_domain=wildcards.target_domain).output.metadata
    output:
        hits = "results/02_sequence_selection/{target_domain}/01_domain_hits/{genome}_protdomain_hits.tsv"
    params:
        e = config["hmmer"]["e"],
        script = "workflow/src/scripts/search_for_domain.py"
    log:
        "logs/02_sequence_selection/{target_domain}/01_domain_hits/{target_domain}_{genome}_protdomain_search.log"
    conda:
        "../envs/domain_analysis_env.yaml"
    threads: 2
    resources:
        mem_mb=config["genomes"]["locate_domains"]["mem_mb"],
        runtime=config["genomes"]["locate_domains"]["runtime"]
    shell:
        """
        hmmsearch --tblout {output.hits} \
                  -E {params.e} \
                  --cpu {threads} \
                  {input.hmm} \
                  {input.protein} \
                  2>&1 | tee {log}
        """

rule extract_flanking_sequences:
    input:
        domain_hit = "results/02_sequence_selection/{target_domain}/01_domain_hits/{genome}_protdomain_hits.tsv",
        genome = "data/genomes/{genome}/genomic.fna",
        annotation = "data/genomes/{genome}/geneannotation.gff"
    output:
        utr5 = "results/02_sequence_selection/{target_domain}/02_flanking_seqs/{genome}_upstream.fasta",
        cds = "results/02_sequence_selection/{target_domain}/02_flanking_seqs/{genome}_cds.fasta",
        utr3 = "results/02_sequence_selection/{target_domain}/02_flanking_seqs/{genome}_downstream.fasta",
        coords = "results/02_sequence_selection/{target_domain}/02_flanking_seqs/{genome}_domain_flanking_coords.bed"
    params:
        script = "workflow/src/extract_flanking_intergenic.py",
        upstream = config["flanking"]["upstream"],
        downstream = config["flanking"]["downstream"],
        temp_dir = ".temp",
        cds_overlap = config["flanking"]["cds_overlap"]
    log:
        "logs/02_sequence_selection/{target_domain}/02_flanking_seqs/{target_domain}_{genome}_flanking_sequence_extraction.log"
    conda:
        "../envs/domain_analysis_env.yaml"
    resources:
        mem_mb=config["genomes"]["extract_flanking_sequences"]["mem_mb"],
        runtime=config["genomes"]["extract_flanking_sequences"]["runtime"]
    shell:
        """
        mkdir -p .temp
        python {params.script} \
            --input-hits {input.domain_hit} \
            --genome {input.genome} \
            --annotation {input.annotation} \
            --output-upstream {output.utr5} \
            --output-cds {output.cds} \
            --output-downstream {output.utr3} \
            --output-coords {output.coords} \
            --upstream {params.upstream} \
            --downstream {params.downstream} \
            --cds-overlap {params.cds_overlap} \
            --temp-dir {params.temp_dir} \
            2>&1 | tee {log}
        """

rule combine_flanking_sequences:
    input:
        utr5 = lambda wildcards: expand(
            "results/02_sequence_selection/{target_domain}/02_flanking_seqs/{accession}_upstream.fasta",
            accession=get_accessions(wildcards),
            target_domain=wildcards.target_domain
        ),
        cds = lambda wildcards: expand(
            "results/02_sequence_selection/{target_domain}/02_flanking_seqs/{accession}_cds.fasta",
            accession=get_accessions(wildcards),
            target_domain=wildcards.target_domain
        ),
        utr3 = lambda wildcards: expand(
            "results/02_sequence_selection/{target_domain}/02_flanking_seqs/{accession}_downstream.fasta",
            accession=get_accessions(wildcards),
            target_domain=wildcards.target_domain
        )
    output:
        utr5 = "results/02_sequence_selection/{target_domain}/03_combined/{target_domain}_domain_flanking_upstream.fasta",
        cds = "results/02_sequence_selection/{target_domain}/03_combined/{target_domain}_domain_flanking_cds.fasta",
        utr3 = "results/02_sequence_selection/{target_domain}/03_combined/{target_domain}_domain_flanking_downstream.fasta"
    log:
        "logs/02_sequence_selection/03_combined/{target_domain}/{target_domain}_flanking_sequence_combination.log"
    conda:
        "../envs/domain_analysis_env.yaml"
    resources:
        mem_mb=config["genomes"]["combine_flanking_sequences"]["mem_mb"],
        runtime=config["genomes"]["combine_flanking_sequences"]["runtime"]
    shell:
        """
        cat {input.utr5} > {output.utr5}
        cat {input.cds} > {output.cds}
        cat {input.utr3} > {output.utr3}
        """

rule filter_flanking_sequences:
    input:
        sequences = "results/02_sequence_selection/{target_domain}/03_combined/{target_domain}_domain_flanking_{region}.fasta"
    output:
        filtered = "results/02_sequence_selection/{target_domain}/03_combined/{target_domain}_domain_flanking_{region}_filtered.fasta"
    params:
        script = "workflow/src/filter_seqs.py",
        length_threshold = config["flanking"]["length_threshold"],
        gc_content_threshold = config["flanking"]["gc_content_threshold"]
    log:
        "logs/02_sequence_selection/03_combined/{target_domain}/filter_flanking_sequences_{region}.log"
    conda:
        "../envs/domain_analysis_env.yaml"
    resources:
        mem_mb=config["genomes"]["filter_flanking_sequences"]["mem_mb"],
        runtime=config["genomes"]["filter_flanking_sequences"]["runtime"]
    shell:
        """
        python {params.script} \
            --input {input.sequences} \
            --output {output.filtered} \
            --length-threshold {params.length_threshold} \
            --gc-threshold {params.gc_content_threshold} \
            2>&1 | tee {log}
        """