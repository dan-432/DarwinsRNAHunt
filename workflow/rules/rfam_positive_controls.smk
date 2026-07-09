


def get_rfam_motif_ids(wildcards):
    return config["controls"]["rfam_families"]

checkpoint get_accession_ids:
    input:
        "resources/rfam/Rfam.full_region"
    output:
        accession_ids = "resources/rfam/{family}_taxonomy.json"
    params:
        script = "workflow/src/rfam_get_accession_ids.py",
        email = config["databases"]["entrez"]["email"]
    conda:
        "../envs/domain_taxonomy_env.yaml"
    log:
        "logs/rfam/{family}/get_accession_ids.log"
    resources:
        mem_mb=8000,  # ~8 GB to safely hold Rfam file + overhead
        disk_mb=5000
    shell:
        """
        python {params.script} \
            --input {input} \
            --family {wildcards.family} \
            --output-list {output.accession_ids} \
            --entrez-email {params.email} \
            > {log} 2>&1
        """

rule select_control_genomes:
    input:
        tree = "resources/gtdb/bac120.tree",
        metadata = "resources/gtdb/bac120_metadata.tsv.gz",
        to_keep = "resources/rfam/{family}_taxonomy.json"
    output:
        tree = "results/00_controls/rfam/{family}/01_phylogeny/tree.nwk",
        metadata = "results/00_controls/rfam/{family}/01_phylogeny/filtered_metadata.tsv.gz",
        genomes = "results/00_controls/rfam/{family}/01_phylogeny/select_genome_accessions.json"
    params:
        script = "workflow/src/extract_subtree.py",
        completeness = config["taxonomy"]["completeness_threshold"],
        contamination = config["taxonomy"]["contamination_threshold"]
    log:
        "logs/rfam/{family}/select_control_genomes.log"
    conda:
        "../envs/domain_taxonomy_env.yaml"
    shell:
        """
        python {params.script} \
            --input-tree {input.tree} \
            --metadata {input.metadata} \
            --to-keep {input.to_keep} \
            --completeness {params.completeness} \
            --contamination {params.contamination} \
            --output-tree {output.tree} \
            --output-metadata {output.metadata} \
            --output-genomes {output.genomes} \
            2>&1 | tee {log}
        """

def get_motif_accession_ids(wildcards):
    from json import load
    # Pass the family wildcard to the checkpoint
    checkpoint_output = checkpoints.download_rfam_motif_genomes.get(family=wildcards.family)
    with open(checkpoint_output.output.metadata) as f:
        metadata = load(f)
    return metadata["accessions"]

checkpoint download_rfam_motif_genomes:
    input:
        genome_accessions = "results/00_controls/rfam/{family}/01_phylogeny/select_genome_accessions.json"
    output:
        metadata = "results/00_controls/rfam/{family}/01_phylogeny/genomes_available_metadata.json"
    params:
        script = "workflow/src/download_ncbi_genomes.sh",
        outdir = "resources/genomes",
        formats = "gff3,genome,protein",
        batch_size = config["databases"]["ncbi"]["batch_size"],
        genome_log = "resources/genomes/genome_log.json",
        temp_dir = ".temp"
    log:
        "logs/genomes/download_genomes_{family}.log"
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

rule extract_rfam_threshold:
    input:
        family_txt = "resources/rfam/database_files/family.txt"
    output:
        threshold = "resources/rfam/cmfind_thresholds/{family}.txt"
    shell:
        """
        THRESHOLD=$(grep "^{wildcards.family}" {input.family_txt} | cut -f 14 | \
                    awk '{{for(i=1;i<=NF;i++) if($i=="-T") print $(i+1)}}')
        echo "$THRESHOLD" > {output.threshold}
        """

rule find_motif:
    input:
        flag = lambda wildcards: checkpoints.download_rfam_motif_genomes.get(family=wildcards.family).output.metadata,
        fasta = "resources/genomes/{genome}/genomic.fna",
        cm_model = "resources/rfam/cm_models/{family}.cm",
        threshold = "resources/rfam/cmfind_thresholds/{family}.txt"
    output:
        "results/00_controls/rfam/{family}/02_flanking/{genome}/motifs.tsv"
    params:
        temp_dir = ".temp"
    threads: config["cmfinder"]["threads"]["initial_discovery"]
    conda:
        "../envs/rna_motif_env.yaml"
    log:
        "logs/rfam/{family}/{genome}/find_motif.log"
    shell:
        """
        set -euo pipefail
        
        mkdir -p $(dirname {output})
        THRESHOLD=$(cat {input.threshold})
        echo "Using threshold: $THRESHOLD" >> {log}
        
        # Work in tmpdir, copy result to NFS
        TMPDIR=$(mktemp -d {params.temp_dir}/find_motif.XXXXXX)
        trap "rm -rf $TMPDIR" EXIT
        
        cmsearch --tblout $TMPDIR/motifs.tsv.raw \
            --cpu {threads} --verbose --nohmmonly \
            -T $THRESHOLD \
            {input.cm_model} {input.fasta} >> {log} 2>&1
        
        awk 'BEGIN {{OFS="\t"}} !/^#/ && NF>=16 {{print $1, $8, $9, $15}}' \
            $TMPDIR/motifs.tsv.raw > $TMPDIR/motifs.tsv 2>> {log}
        
        # Copy final result to output (bypasses NFS touch issues)
        mv $TMPDIR/motifs.tsv {output}
        """

rule find_flanking_cds:
    input:
        motif_index = "results/00_controls/rfam/{family}/02_flanking/{genome}/motifs.tsv",
        gff = "resources/genomes/{genome}/geneannotation.gff",
        proteome = "resources/genomes/{genome}/protein.faa"
    output:
        "results/00_controls/rfam/{family}/02_flanking/{genome}/flanking_cds.faa"
    params:
        script = "workflow/src/extract_flanking_cds.py",
        search_area = config["controls"]["cds_search_area"],
        genome_accession = "{genome}"
    conda:
        "../envs/domain_analysis_env.yaml"
    log:
        "logs/rfam/{family}/{genome}/find_flanking_cds.log"
    shell:
        """
        python {params.script} \
            --regions {input.motif_index} \
            --family {wildcards.family} \
            --gff {input.gff} \
            --net-size {params.search_area} \
            --assembly-accession {params.genome_accession} \
            --proteome {input.proteome} \
            --output {output} \
            2>&1 | tee {log}
        """

rule combine_flanking_cds:
    input:
        lambda wildcards: expand(
            "results/00_controls/rfam/{family}/02_flanking/{genome}/flanking_cds.faa",
            genome=get_motif_accession_ids(wildcards),
            family=wildcards.family
        )
    output:
        "results/00_controls/rfam/{family}/02_flanking/combined_flanking_cds.faa"
    run:
        with open(output[0], "w") as outfile:
            for infile in input:
                with open(infile) as f:
                    outfile.write(f.read())

rule create_mmseq2_db:
    input:
        "results/00_controls/rfam/{family}/02_flanking/combined_flanking_cds.faa"
    output:
        directory("results/00_controls/rfam/{family}/03_mmseq_dbs/combined_flanking_cds_db")
    params:
        db_prefix = "db"
    conda:
        "../envs/rna_motif_env.yaml"
    log:
        "logs/rfam/{family}/cluster_flanking_cds_db.log"
    shell:
        """
        mkdir -p {output}
        mmseqs createdb {input} {output}/{params.db_prefix} 2>&1 | tee {log}
        """

def get_cluster_ids(wildcards):
    """
    Aggregate cluster IDs from checkpoint output.
    Only return clusters with >1 sequence (skip singletons).
    """
    clusters_file = checkpoints.cluster_flanking_cds.get(**wildcards).output.table
    
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

checkpoint cluster_flanking_cds:
    input:
        "results/00_controls/rfam/{family}/03_mmseq_dbs/combined_flanking_cds_db"
    output:
        table = "results/00_controls/rfam/{family}/04_clusters/flanking_cds_clusters.tsv",
        cluster_db = directory("results/00_controls/rfam/{family}/03_mmseq_dbs/flanking_cds_clustering_db")
    params:
        identity_threshold = config["clustering"]["mmseqs2"]["min_seq_id"],
        coverage_threshold = config["clustering"]["mmseqs2"]["coverage_threshold"],
        cov_mode = config["clustering"]["mmseqs2"]["cov_mode"],
        cluster_mode = config["clustering"]["mmseqs2"]["cluster_mode"],
        sensitivity = config["clustering"]["mmseqs2"]["sensitivity"],
        temp_dir = ".temp",
        cluster_db_prefix = "cluster_db",
        flanking_db_prefix = "db",
        refactor_script = "workflow/src/remap_cluster_ids.py"
    threads:
        config["clustering"]["mmseqs2"]["threads"]
    conda:
        "../envs/rna_motif_env.yaml"
    log:
        "logs/rfam/{family}/cluster_flanking_cds.log"
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

# rule get_family_clustering:
#     input:
#         clusters = lambda wildcards: expand(
#             "results/00_controls/rfam/{family}/flanking_cds_clusters.tsv",
#             family=get_rfam_motif_ids(wildcards)
#         )

#     output:
#         "results/00_controls/rfam/family_clustering_summary.tsv"
#     run:
#         import pandas as pd
#         cluster_dfs = []
#         for infile in input:
#             df = pd.read_csv(infile, sep="\t", header=None, names=["ClusterID", "SequenceID"])
#             df["Family"] = infile.split("/")[3]  # Extract family from path
#             cluster_dfs.append(df)
#         combined_df = pd.concat(cluster_dfs)
#         summary = combined_df.groupby("Family").agg({
#             "ClusterID": "nunique",
#             "SequenceID": "nunique"
#         }).reset_index().rename(columns={"ClusterID": "NumClusters", "SequenceID": "NumSequences"})
#         summary.to_csv(output[0], sep="\t", index=False)


rule extract_cluster_fasta:
    """Extract protein sequences for each cluster."""
    input:
        fasta="results/00_controls/rfam/{family}/02_flanking/combined_flanking_cds.faa",
        clusters="results/00_controls/rfam/{family}/04_clusters/flanking_cds_clusters.tsv"
    output:
        cluster_fasta="results/00_controls/rfam/{family}/04_clusters/cluster_{cluster_id}.faa"
    params:
        script="workflow/src/extract_cluster_fasta.py"
    conda:
        "../envs/domain_analysis_env.yaml"
    log:
        "logs/rfam/{family}/extract_cluster_{cluster_id}.log"
    shell:
        """
        python {params.script} \
            --fasta {input.fasta} \
            --clusters {input.clusters} \
            --cluster_id '{wildcards.cluster_id}' \
            --output {output.cluster_fasta} \
            2>&1 | tee {log}
        """

rule mafft_align:
    """Align cluster proteins using MAFFT."""
    input:
        fasta="results/00_controls/rfam/{family}/04_clusters/cluster_{cluster_id}.faa"
    output:
        msa="results/00_controls/rfam/{family}/05_alignments/cluster_{cluster_id}.msa"
    conda:
        "../envs/rna_motif_env.yaml"
    threads: 4
    log:
        "logs/rfam/{family}/mafft_cluster_{cluster_id}.log"
    shell:
        """
        mafft --thread {threads} --auto {input.fasta} > {output.msa} 2> {log}
        """

rule hmmer_build:
    """Build HMM profile from MSA using HMMER."""
    input:
        msa=ancient("results/00_controls/rfam/{family}/05_alignments/cluster_{cluster_id}.msa")
    output:
        hmm="results/00_controls/rfam/{family}/06_hmm_profiles/cluster_{cluster_id}.hmm"
    conda:
        "../envs/domain_analysis_env.yaml"
    log:
        "logs/rfam/{family}/hmmbuild_cluster_{cluster_id}.log"
    shell:
        """
        hmmbuild --amino {output.hmm} {input.msa} 2>&1 | tee {log}
        sleep 10
        touch {output.hmm}
        """

rule hmmcalibrate:
    """Calibrate HMM for statistical significance testing."""
    input:
        hmm="results/00_controls/rfam/{family}/06_hmm_profiles/cluster_{cluster_id}.hmm"
    output:
        hmm_calib="results/00_controls/rfam/{family}/06_hmm_profiles/cluster_{cluster_id}.hmm"
    conda:
        "../envs/domain_analysis_env.yaml"
    log:
        "logs/rfam/{family}/hmmcalibrate_cluster_{cluster_id}.log"
    shell:
        """
        hmmcalibrate {input.hmm} 2>&1 | tee {log}
        """

rule cluster_spread_analysis:
    input:
        clusters = "results/00_controls/rfam/{family}/04_clusters/flanking_cds_clusters.tsv",
        fasta = "results/00_controls/rfam/{family}/02_flanking/combined_flanking_cds.faa",
        metadata = "results/00_controls/rfam/{family}/01_phylogeny/genomes_available_metadata.json",
        hmms = lambda wildcards: expand(
            "results/00_controls/rfam/{family}/06_hmm_profiles/cluster_{cluster_id}.hmm",
            cluster_id=get_cluster_ids(wildcards), family=wildcards.family
        )
    output:
        report = "results/00_controls/rfam/{family}/07_cluster_analysis/spread_report.tsv",
        presence_matrix = "results/00_controls/rfam/{family}/07_cluster_analysis/presence_absence_matrix.tsv",
        summary = "results/00_controls/rfam/{family}/07_cluster_analysis/summary.txt"
    params:
        script = "workflow/src/analyse_cluster_spread.py"
    conda:
        "../envs/domain_analysis_env.yaml"
    log:
        "logs/rfam/{family}/cluster_spread_analysis.log"
    shell:
        """
        python {params.script} \
            --clusters {input.clusters} \
            --fasta {input.fasta} \
            --metadata {input.metadata} \
            --output_report {output.report} \
            --output_presence_matrix {output.presence_matrix} \
            --output_summary {output.summary} \
            2>&1 | tee {log}
        """

rule analyse_cluster_phylogeny:
    input:
        presence_matrix = "results/00_controls/rfam/{family}/07_cluster_analysis/presence_absence_matrix.tsv",
        tree = "results/00_controls/rfam/{family}/01_phylogeny/tree.nwk",
        metadata = "results/00_controls/rfam/{family}/01_phylogeny/filtered_metadata.tsv.gz"
    output:
        report = "results/00_controls/rfam/{family}/07_cluster_analysis/phylo_report.tsv",
        itol = "results/00_controls/rfam/{family}/07_cluster_analysis/itol_presence.txt"
    params:
        script = "workflow/src/analyse_cluster_phylogeny.py",
        top_n = 20,
        max_unmatched = 300
    conda:
        "../envs/domain_taxonomy_env.yaml"
    log:
        "logs/rfam/{family}/analyse_cluster_phylogeny.log"
    shell:
        """
        python {params.script} \
            --presence-matrix {input.presence_matrix} \
            --tree {input.tree} \
            --tree-metadata {input.metadata} \
            --output-report {output.report} \
            --output-itol {output.itol} \
            --top-n {params.top_n} \
            --expected-max-unmatched {params.max_unmatched} \
            2>&1 | tee {log}
        """

rule get_family_clustering:
    input:
        clusters = lambda wildcards: expand(
            "results/00_controls/rfam/{family}/07_cluster_analysis/phylo_report.tsv",
            family=get_rfam_motif_ids(wildcards)
        )

    output:
        "results/00_controls/rfam/family_phylo_summary.tsv"
    shell:
        """
        cat {input.clusters} > {output}"""