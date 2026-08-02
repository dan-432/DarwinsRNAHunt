""" 
Pipeline for motif discovery by multiple sequence allignment and R-scape.

first cluster with mmseqs2 then align clusters with mafft then check for conservd structure in allignment with R-scape

"""

rule create_motif_mmseq2_db:
    input:
        "results/02_sequence_selection/{target_domain}/03_combined/{target_domain}_domain_flanking_{region}_filtered.fasta"
    output:
        directory("results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/01_clustering/combined_flanking_dbs")
    params:
        db_prefix = "db"
    conda:
        "../envs/rna_motif_env.yaml"
    log:
        "logs/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/01_clustering/cluster_flanking_cds_db.log"
    shell:
        """
        mkdir -p {output}
        mmseqs createdb {input} {output}/{params.db_prefix} 2>&1 | tee {log}
        """

def get_motif_cluster_ids(wildcards):
    """
    Aggregate cluster IDs from checkpoint output.
    Only return clusters with >1 sequence (skip singletons).
    """
    clusters_file = checkpoints.cluster_flanking_regions.get(**wildcards).output.table
    
    cluster_counts = {}
    try:
        with open(clusters_file) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    cluster_id = parts[0]  # Column 1 is the cluster ID
                    cluster_counts[cluster_id] = cluster_counts.get(cluster_id, 0) + 1
    except FileNotFoundError:
        return []
    
    # Only return clusters with >1 sequence
    return sorted([cid for cid, count in cluster_counts.items() if count > config["controls"]["min_cluster_size"]])

checkpoint cluster_flanking_regions:
    input:
        "results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/01_clustering/combined_flanking_dbs"
    output:
        table = "results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/01_clustering/flanking_cds_clusters.tsv",
        cluster_db = directory("results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/01_clustering/clustering_db")
    params:
        identity_threshold = config["motif_clustering"]["mmseqs2"]["min_seq_id"],
        coverage_threshold = config["motif_clustering"]["mmseqs2"]["coverage_threshold"],
        cov_mode = config["motif_clustering"]["mmseqs2"]["cov_mode"],
        cluster_mode = config["motif_clustering"]["mmseqs2"]["cluster_mode"],
        sensitivity = config["motif_clustering"]["mmseqs2"]["sensitivity"],
        temp_dir = ".temp",
        cluster_db_prefix = "cluster_db",
        flanking_db_prefix = "db",
        refactor_script = "workflow/src/remap_cluster_ids.py"
    threads:
        config["motif_clustering"]["mmseqs2"]["threads"]
    conda:
        "../envs/rna_motif_env.yaml"
    log:
        "logs/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/01_clustering/cluster_flanking_cds.log"
    shell:
        """
        mkdir -p {params.temp_dir}
        mkdir -p {output.cluster_db}
        
        mmseqs cluster {input}/{params.flanking_db_prefix} {output.cluster_db}/{params.cluster_db_prefix} {params.temp_dir} \
            --min-seq-id {params.identity_threshold} -c {params.coverage_threshold} --cov-mode {params.cov_mode} \
            --cluster-mode {params.cluster_mode} -s {params.sensitivity} 2>&1 | tee {log}
        
        mmseqs createtsv {input}/{params.flanking_db_prefix} {input}/{params.flanking_db_prefix} {output.cluster_db}/{params.cluster_db_prefix} {output.table}.raw 2>&1 | tee -a {log}
        
        python {params.refactor_script} --input {output.table}.raw --output {output.table} 2>&1 | tee -a {log}
        """

rule extract_motif_cand_cluster_fasta:
    """Extract protein sequences for each cluster."""
    input:
        fasta="results/02_sequence_selection/{target_domain}/03_combined/{target_domain}_domain_flanking_{region}_filtered.fasta",
        clusters="results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/01_clustering/flanking_cds_clusters.tsv"
    output:
        cluster_fasta="results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/01_clustering/cluster_{cluster_id}.fasta"
    params:
        script="workflow/src/extract_cluster_nucleotide_fasta.py"
    conda:
        "../envs/domain_analysis_env.yaml"
    log:
        "logs/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/01_clustering/extract_cluster_{cluster_id}.log"
    shell:
        """
        python {params.script} \
            --fasta {input.fasta} \
            --clusters {input.clusters} \
            --cluster_id '{wildcards.cluster_id}' \
            --output {output.cluster_fasta} \
            2>&1 | tee {log}
        """

rule mafft_align_motif_cand:
    """Align cluster proteins using MAFFT."""
    input:
        fasta="results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/01_clustering/cluster_{cluster_id}.fasta"
    output:
        msa="results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/02_alignments/cluster_{cluster_id}.msa"
    conda:
        "../envs/rna_motif_env.yaml"
    threads: 4
    log:
        "logs/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/02_alignments/mafft_cluster_{cluster_id}.log"
    shell:
        """
        mafft --thread {threads} --auto {input.fasta} > {output.msa} 2> {log}
        """

rule build_stockholm:
    """convert msa output to stockholm format with esl"""
    input:
        msa="results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/02_alignments/cluster_{cluster_id}.msa"
    output:
        sto="results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/02_alignments/cluster_{cluster_id}.sto"
    conda:
        "../envs/rna_motif_env.yaml"
    threads: 1
    log:
        "logs/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/02_alignments/esl_reformat_{cluster_id}.log"
    shell:
        """
        esl-reformat stockholm {input.msa} > {output.sto} 2> {log}
        """

rule find_structure:
    """Score conserved structure with R-scape"""
    input:
        sto="results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/02_alignments/cluster_{cluster_id}.sto"
    output:
        struc_tbl = "results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/03_structure/cluster_{cluster_id}_struc.txt"
    conda:
        "../envs/rna_motif_env.yaml"
    log:
        "logs/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/03_structure/cluster_{cluster_id}_struc.log"
    shell:
        """
        rscape --cacofold  {input.sto} > {output.struc_tbl} 2> {log}
        """


rule aggregate_structure_results:
    """Resolve clusters, can later implement some sort of combined report"""
    input:
        struc_tbls = lambda wildcards: expand( 
            "results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/03_structure/cluster_{cluster_id}_struc.txt", 
            target_domain = wildcards.target_domain, region = wildcards.region,
            cluster_id=get_motif_cluster_ids(wildcards))
    output:
        report = "results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/final_report.txt"
    log:
        "logs/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/final_report.log"
    shell:
        "cat {input.struc_tbls} > {output.report} 2> {log}"