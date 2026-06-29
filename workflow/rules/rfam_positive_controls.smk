


def get_rfam_motif_ids(wildcards):
    import json
    with open("results/00_controls/rfam/select_motif_controls.json") as f:
        data = json.load(f)
    return data

checkpoint get_accession_ids:
    input:
        "resources/rfam/Rfam.full_region"
    output:
        accession_ids = "results/00_controls/rfam/{family}/accession_ids.json"
    params:
        script = "workflow/src/rfam_get_accession_ids.py",
        email = config["databases"]["entrez"]["email"]
    conda:
        "../envs/domain_taxonomy_env.yaml"
    log:
        "logs/rfam/{family}/get_accession_ids.log"
    shell:
        """
        python {params.script} \
            --input {input} \
            --family {wildcards.family} \
            --output-list {output.accession_ids} \
            --entrez-email {params.email} \
            > {log} 2>&1
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
        genome_accessions = "results/00_controls/rfam/{family}/accession_ids.json"
    output:
        metadata = "resources/rfam_{family}_genomes/download_metadata.json",
        flag = "resources/rfam_{family}_genomes/.ncbi_download_complete.flag"
    params:
        script = "workflow/src/download_ncbi_genomes.sh",
        outdir = "resources/rfam_{family}_genomes",
        formats = "gff3,genome,protein",
        batch_size = config["databases"]["ncbi"]["batch_size"]
    log:
        "logs/genomes/download_genomes_{family}.log"
    conda:
        "../envs/genome_download_env.yaml"
    shell:
        """
        bash {params.script} \
            {input.genome_accessions} \
            {params.outdir} \
            {params.formats} \
            {params.batch_size} \
            2>&1 | tee {log}
        
        touch {output.flag}
        """

rule find_motif:
    input:
        flag = "resources/rfam_{family}_genomes/.ncbi_download_complete.flag",
        fasta = "resources/rfam_{family}_genomes/{genome}/genomic.fna",
        family_db = "resources/rfam/database_files/family.txt",
        cm_model = "resources/rfam/cm_models/{family}.cm"
    output:
        "results/00_controls/rfam/{family}/flanking/{genome}/motifs.tsv"
    threads: config["cmfinder"]["threads"]["initial_discovery"]
    conda:
        "../envs/rna_motif_env.yaml"
    log:
        "logs/rfam/{family}/{genome}/find_motif.log",
    shell:
        """
        # Extract -T threshold from family.txt
        THRESHOLD=$(grep "^{wildcards.family}" {input.family_db} | cut -f 14 | \
            awk '{{for(i=1;i<=NF;i++) if($i=="-T") print $(i+1)}}')
        
        echo "Using threshold: $THRESHOLD" >> {log}
        
        # Run cmsearch with extracted threshold
        cmsearch --tblout {output}.raw \
            --cpu {threads} --verbose --nohmmonly \
            -T $THRESHOLD \
            {input.cm_model} {input.fasta} >> {log} 2>&1
        
        # Parse tblout to TSV (seq_id, start, end, bit_score)
        awk 'BEGIN {{OFS="\t"}} 
             !/^#/ && NF>=16 {{print $1, $8, $9, $15}}' \
            {output}.raw > {output}
        
        rm {output}.raw
        """

rule find_flanking_cds:
    input:
        motif_index = "results/00_controls/rfam/{family}/flanking/{genome}/motifs.tsv",
        gff = "resources/rfam_{family}_genomes/{genome}/geneannotation.gff",
        proteome = "resources/rfam_{family}_genomes/{genome}/protein.faa"
    output:
        "results/00_controls/rfam/{family}/flanking/{genome}/flanking_cds.faa"
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
            "results/00_controls/rfam/{family}/flanking/{genome}/flanking_cds.faa",
            genome=get_motif_accession_ids(wildcards),
            family=wildcards.family
        )
    output:
        "results/00_controls/rfam/{family}/combined_flanking_cds.faa"
    run:
        with open(output[0], "w") as outfile:
            for infile in input:
                with open(infile) as f:
                    outfile.write(f.read())

rule create_mmseq2_db:
    input:
        "results/00_controls/rfam/{family}/combined_flanking_cds.faa"
    output:
        "results/00_controls/rfam/{family}/combined_flanking_cds_db"
    conda:
        "../envs/rna_motif_env.yaml"
    log:
        "logs/rfam/{family}/cluster_flanking_cds_db.log"
    shell:
        """
        mmseqs createdb {input} {output}
            2>&1 | tee {log}
        """

def get_cluster_ids(wildcards):
    """
    Aggregate cluster IDs from checkpoint output.
    Called after mmseqs2 clustering to determine how many clusters exist.
    """
    checkpoint_output = checkpoints.mmseqs2_cluster.get(**wildcards).params.cluster_db
    clusters_file = checkpoints.mmseqs2_cluster.get(**wildcards).output.table
    
    cluster_ids = set()
    try:
        with open(clusters_file) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 1:
                    cluster_id = parts[0]
                    cluster_ids.add(cluster_id)
    except FileNotFoundError:
        # Checkpoint hasn't run yet
        return []
    
    return sorted(cluster_ids)

checkpoint cluster_flanking_cds:
    input:
        "results/00_controls/rfam/{family}/combined_flanking_cds_db"
    output:
        table = "results/00_controls/rfam/{family}/flanking_cds_clusters.tsv"
    params:
        identity_threshold = config["clustering"]["mmseqs2"]["min_seq_id"],
        temp_dir = ".temp",
        cluster_db = "results/00_controls/rfam/{family}/flanking_cds_clustering_db"
    conda:
        "../envs/rna_motif_env.yaml"
    log:
        "logs/rfam/{family}/cluster_flanking_cds.log"
    shell:
        """
        mmseqs cluster {input} {params.cluster_db} {params.temp_dir} --min-seq-id {params.identity_threshold} \
            2>&1 | tee {log}
        mmseqs createtsv {input} {input} {params.cluster_db} {output}
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
        fasta="results/00_controls/rfam/{family}/combined_flanking_cds.faa",
        clusters="results/00_controls/rfam/{family}/flanking_cds_clusters.tsv"
    output:
        cluster_fasta="results/00_controls/rfam/{family}/clusters/cluster_{cluster_id}.faa"
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
            --cluster_id {wildcards.cluster_id} \
            --output {output.cluster_fasta} \
            2>&1 | tee {log}
        """

rule mafft_align:
    """Align cluster proteins using MAFFT."""
    input:
        fasta="results/00_controls/rfam/{family}/clusters/cluster_{cluster_id}.faa"
    output:
        msa="results/00_controls/rfam/{family}/alignments/cluster_{cluster_id}.msa"
    conda:
        "../envs/alignment_env.yaml"
    threads: 4
    log:
        "logs/rfam/{family}/mafft_cluster_{cluster_id}.log"
    shell:
        """
        mafft --thread {threads} --auto {input.fasta} > {output.msa} 2>&1 | tee {log}
        """

rule hmmer_build:
    """Build HMM profile from MSA using HMMER."""
    input:
        msa="results/00_controls/rfam/{family}/alignments/cluster_{cluster_id}.msa"
    output:
        hmm="results/00_controls/rfam/{family}/hmm_profiles/cluster_{cluster_id}.hmm"
    conda:
        "../envs/hmmer_env.yaml"
    log:
        "logs/rfam/{family}/hmmbuild_cluster_{cluster_id}.log"
    shell:
        """
        hmmbuild --amino {output.hmm} {input.msa} 2>&1 | tee {log}
        """

rule hmmcalibrate:
    """Calibrate HMM for statistical significance testing."""
    input:
        hmm="results/00_controls/rfam/{family}/hmm_profiles/cluster_{cluster_id}.hmm"
    output:
        hmm_calib="results/00_controls/rfam/{family}/hmm_profiles/cluster_{cluster_id}.hmm"
    conda:
        "../envs/hmmer_env.yaml"
    log:
        "logs/rfam/{family}/hmmcalibrate_cluster_{cluster_id}.log"
    shell:
        """
        hmmcalibrate {input.hmm} 2>&1 | tee {log}
        """

rule cluster_spread_analysis:
    """Analyze how clusters spread across genomes."""
    input:
        clusters="results/00_controls/rfam/{family}/mmseqs2/clusters.tsv",
        fasta="results/00_controls/rfam/{family}/combined_flanking_cds.faa",
        hmms=expand("results/00_controls/rfam/{{family}}/hmm_profiles/cluster_{cluster_id}.hmm",
                    cluster_id=get_cluster_ids)
    output:
        report="results/00_controls/rfam/{family}/cluster_analysis/spread_report.tsv",
        summary="results/00_controls/rfam/{family}/cluster_analysis/summary.txt"
    params:
        script="workflow/src/analyze_cluster_spread.py"
    conda:
        "../envs/domain_analysis_env.yaml"
    log:
        "logs/rfam/{family}/cluster_spread_analysis.log"
    shell:
        """
        python {params.script} \
            --clusters {input.clusters} \
            --fasta {input.fasta} \
            --output_report {output.report} \
            --output_summary {output.summary} \
            2>&1 | tee {log}
        """