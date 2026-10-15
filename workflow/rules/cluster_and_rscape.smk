""" 
Pipeline for motif discovery by multiple sequence allignment and R-scape.

first cluster with mmseqs2 then align clusters with mafft then check for conservd structure in allignment with R-scape

mmseqs2 search (kept, with coordinates) -> mmseqs2 clust (membership only)
-> trim each cluster to its core homologous window -> mafft -> R-scape
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

rule mmseqs2_search:
    """
    All-vs-all local alignment search. This does what the search half of
    `mmseqs cluster` does internally, except we keep the alignment DB so
    coordinates are available for trimming afterward instead of being
    discarded.
    """
    input:
        db_dir = "results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/01_clustering/combined_flanking_dbs"
    output:
        search_db = directory("results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/01_clustering/search_db"),
        alignments = "results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/01_clustering/flanking_cds_alignments.m8"
    params:
        identity_threshold = config["motif_clustering"]["mmseqs2"]["min_seq_id"],
        coverage_threshold = config["motif_clustering"]["mmseqs2"]["coverage_threshold"],
        cov_mode = config["motif_clustering"]["mmseqs2"]["cov_mode"],
        temp_dir = ".temp_search_{target_domain}_{region}",
        search_db_prefix = "search",
        flanking_db_prefix = "db",
        format_output = "query,target,qstart,qend,tstart,tend,qlen,tlen,fident,bits,evalue"
    threads:
        config["motif_clustering"]["threads"]
    conda:
        "../envs/rna_motif_env.yaml"
    log:
        "logs/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/01_clustering/mmseqs_search.log"
    shell:
        """
        mkdir -p {params.temp_dir}
        mkdir -p {output.search_db}

        mmseqs search {input.db_dir}/{params.flanking_db_prefix} {input.db_dir}/{params.flanking_db_prefix} \
            {output.search_db}/{params.search_db_prefix} {params.temp_dir} \
            --min-seq-id {params.identity_threshold} -c {params.coverage_threshold} --cov-mode {params.cov_mode} \
            --search-type 3 \
            -a --alignment-mode 3 --threads {threads} 2>&1 | tee {log}

        mmseqs convertalis {input.db_dir}/{params.flanking_db_prefix} {input.db_dir}/{params.flanking_db_prefix} \
            {output.search_db}/{params.search_db_prefix} {output.alignments} \
            --format-output "{params.format_output}" 2>&1 | tee -a {log}

        rm -rf {params.temp_dir}
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
                    cluster_id = parts[0]
                    cluster_counts[cluster_id] = cluster_counts.get(cluster_id, 0) + 1
    except FileNotFoundError:
        return []

    return sorted([cid for cid, count in cluster_counts.items() if count > config["motif_clustering"]["min_cluster_size"]])

checkpoint cluster_flanking_regions:
    """
    Cluster using the already-computed search results, rather than
    re-running search internally as `mmseqs cluster` would.
    """
    input:
        db_dir = "results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/01_clustering/combined_flanking_dbs",
        search_db = "results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/01_clustering/search_db"
    output:
        table = "results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/01_clustering/flanking_cds_clusters.tsv",
        cluster_db = directory("results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/01_clustering/clustering_db")
    params:
        cluster_mode = config["motif_clustering"]["mmseqs2"]["cluster_mode"],
        cluster_db_prefix = "cluster_db",
        search_db_prefix = "search",
        flanking_db_prefix = "db",
        refactor_script = "workflow/src/remap_cluster_ids.py"
    conda:
        "../envs/rna_motif_env.yaml"
    log:
        "logs/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/01_clustering/cluster_flanking_cds.log"
    shell:
        """
        mkdir -p {output.cluster_db}

        mmseqs clust {input.db_dir}/{params.flanking_db_prefix} {input.search_db}/{params.search_db_prefix} \
            {output.cluster_db}/{params.cluster_db_prefix} --cluster-mode {params.cluster_mode} 2>&1 | tee {log}

        mmseqs createtsv {input.db_dir}/{params.flanking_db_prefix} {input.db_dir}/{params.flanking_db_prefix} \
            {output.cluster_db}/{params.cluster_db_prefix} {output.table}.raw 2>&1 | tee -a {log}

        python {params.refactor_script} --input {output.table}.raw --output {output.table} 2>&1 | tee -a {log}
        """

rule extract_motif_cand_cluster_fasta:
    """Extract nucleotide sequences for each cluster (unchanged)."""
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

rule trim_cluster_to_core_window:
    """
    Trim every member of a cluster down to the region that locally aligns
    to an in-cluster anchor sequence, using coordinates already computed
    by the single mmseqs2_search pass (no per-cluster re-search needed).
    """
    input:
        fasta="results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/01_clustering/cluster_{cluster_id}.fasta",
        alignments="results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/01_clustering/flanking_cds_alignments.m8"
    output:
        trimmed="results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/01_clustering/trimmed_cluster_{cluster_id}.fasta"
    params:
        script="workflow/src/trim_to_conserved_window.py",
        padding=config["trimming"]["padding_nt"],
        min_len=config["trimming"]["min_trim_len"],
        max_len=config["trimming"]["max_trim_len"]
    conda:
        "../envs/domain_analysis_env.yaml"
    log:
        "logs/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/01_clustering/trim_{cluster_id}.log"
    shell:
        """
        python {params.script} \
            --fasta {input.fasta} \
            --alignments {input.alignments} \
            --padding {params.padding} \
            --min-length {params.min_len} \
            --max-length {params.max_len} \
            --output {output.trimmed} \
            2>&1 | tee {log}
        """

rule mafft_align_motif_cand:
    """Align the trimmed, comparably-sized cluster members."""
    input:
        fasta="results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/01_clustering/trimmed_cluster_{cluster_id}.fasta"
    output:
        msa="results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/02_alignments/cluster_{cluster_id}.msa"
    conda:
        "../envs/rna_motif_env.yaml"
    threads: config["motif_clustering"]["mafft"]["threads"]
    log:
        "logs/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/02_alignments/mafft_cluster_{cluster_id}.log"
    shell:
        """
        mafft --thread {threads} --localpair --maxiterate 1000 {input.fasta} > {output.msa} 2> {log}
        """

# build_stockholm, find_structure, build_cm_from_msa, calibrate_cm_model,
# aggregate_structure_results are unchanged from your current file.

rule build_stockholm:
    """convert msa output to stockholm format with esl"""
    input:
        msa="results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/02_alignments/cluster_{cluster_id}.msa"
    output:
        sto="results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/02_alignments/cluster_{cluster_id}.sto"
    conda:
        "../envs/rna_motif_env.yaml"
    threads: config["motif_clustering"]["esl"]["threads"]
    log:
        "logs/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/02_alignments/esl_reformat_{cluster_id}.log"
    shell:
        """
        if [ ! -s {input.msa} ]; then
            echo "Error: Input MSA file {input.msa} is empty or does not exist." >> {log}
            touch {output.sto}  # Create an empty Stockholm file to avoid downstream errors
            exit 0
        fi
        esl-reformat stockholm {input.msa} > {output.sto} 2> {log}
        """

rule find_structure:
    """Score conserved structure with R-scape"""
    input:
        sto="results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/02_alignments/cluster_{cluster_id}.sto"
    output:
        struc_tbl = "results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/03_structure/cluster_{cluster_id}_struc.txt",
        output_dir = directory("results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/03_structure/cluster_{cluster_id}_rscape_output")
    conda:
        "../envs/rna_motif_env.yaml"
    log:
        "logs/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/03_structure/cluster_{cluster_id}_struc.log"
    shell:
        """
        mkdir -p {output.output_dir}
        R-scape --cacofold --outdir {output.output_dir} {input.sto} > {output.struc_tbl} 2> {log} \
            || {{ echo "R-scape failed on cluster {wildcards.cluster_id}, likely too little data for null fit" >> {log}; touch {output.struc_tbl}; }}
        """

def has_covarying_pair(cov_path):
    """True if an R-scape .cov file reports at least one significant covarying pair."""
    if not os.path.exists(cov_path) or os.path.getsize(cov_path) == 0:
        return False
    with open(cov_path) as f:
        lines = [l for l in f if l.strip() and not l.startswith("#")]
    if any("no significant pairs" in l for l in lines):
        return False
    return len(lines) > 0

def get_covarying_cluster_ids(wildcards):
    table = checkpoints.filter_covarying_clusters.get(
        target_domain=wildcards.target_domain, region=wildcards.region
    ).output.table
    with open(table) as f:
        return [line.strip() for line in f if line.strip()]

checkpoint filter_covarying_clusters:
    """
    Runs find_structure for every candidate cluster (ordinary expand
    dependency, no checkpoint needed there), then records which ones
    actually show significant covariation.
    """
    input:
        rscape_dirs = lambda wildcards: expand(
            "results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/03_structure/cluster_{cluster_id}_rscape_output",
            target_domain=wildcards.target_domain, region=wildcards.region,
            cluster_id=get_motif_cluster_ids(wildcards))
    output:
        table = "results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/03_structure/covarying_clusters.tsv"
    run:
        cluster_ids = get_motif_cluster_ids(wildcards)
        with open(output.table, "w") as out:
            for cluster_id, rscape_dir in zip(cluster_ids, input.rscape_dirs):
                cov_file = os.path.join(rscape_dir, f"cluster_{cluster_id}_1.cacofold.cov")
                if has_covarying_pair(cov_file):
                    out.write(f"{cluster_id}\n")

rule build_cm_from_msa:
    input:
        sto = "results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/03_structure/cluster_{cluster_id}_rscape_output"
    output:
        cm = "results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/04_cms/{cluster_id}_uncalibrated.cm",
    conda:
        "../envs/rna_motif_env.yaml"
    log:
        "logs/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/03_structure/cmbuild_{cluster_id}.log"
    shell:
        """
        if [ ! -f "{input.sto}/cluster_{wildcards.cluster_id}_1.cacofold.sto" ]; then
            echo "Error: Stockholm file {input.sto}/cluster_{wildcards.cluster_id}_1.cacofold.sto not found." >> {log}
            touch {output.cm}  # Create an empty CM file to avoid downstream errors
            exit 0
        fi
        echo "Building CM with cmbuild..." >> {log}
        cmbuild -F {output.cm} {input.sto}/cluster_{wildcards.cluster_id}_1.cacofold.sto 2>&1 | tee -a {log}
        echo "CM building complete for {wildcards.cluster_id}" >> {log} 
        """

# Calibrate covariance models (OPTIONAL - can skip for faster results)
rule calibrate_cm_model:
    input:
        cm = "results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/04_cms/{cluster_id}_uncalibrated.cm",
        rscape_out = "results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/03_structure/cluster_{cluster_id}_rscape_output"
    output:
        cm = "results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/04_cms/{cluster_id}_calibrated.cm",
        cm_stats = "results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/04_cms/{cluster_id}_calibrated.cm.stats"
    params:
        tailp = config["cmcalibrate"]["tailp"],
        max_est_runtime = config["infernal"]["cmcalibrate"]["max_runtime"],
        total_length_random_seq = config["cmcalibrate"]["total_length_random_seq"],
        model_name = "cluster_{cluster_id}_1.cacofold"
    conda:
        "../envs/rna_motif_env.yaml"
    threads: config["infernal"]["cmcalibrate"]["threads"]
    resources:
        mem_mb = config["infernal"]["cmcalibrate"]["mem_mb"],
        mem_safety_factor = config["infernal"]["cmcalibrate"]["mem_safety_factor"]
    log:
        "logs/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/04_cms/cmcalibrate_{cluster_id}.log"
    shell:
        """
        set -euo pipefail
        echo "Calibrating motif {wildcards.cluster_id}" > {log}
        if [[ ! -s {input.cm} ]] || [[ ! -s {input.rscape_out}/cluster_{wildcards.cluster_id}_1.cacofold.cov ]] || grep -q "no significant pairs" {input.rscape_out}/cluster_{wildcards.cluster_id}_1.cacofold.cov; then
            echo "Error: Input CM file {input.cm} is empty or does not exist." >> {log}
            touch {output.cm}  # Create an empty CM file to avoid downstream errors
            touch {output.cm_stats}  # Create an empty CM stats file to avoid downstream errors
            exit 0
        fi

        # Copy CM to calibration directory
        cp {input.cm} {output.cm}.tmp

        estimate_output=$(mktemp)
        trap 'rm -f "$estimate_output"' EXIT
        if ! cmcalibrate --memreq --forecast --nforecast {threads} \
            -L {params.total_length_random_seq} {output.cm}.tmp \
            > "$estimate_output" 2>&1; then
            echo "cmcalibrate estimate failed; leaving model uncalibrated." >> {log}
            cat "$estimate_output" >> {log}
            cmstat {input.cm} >> {output.cm_stats} 2>&1 | tee -a {log}
            mv {output.cm}.tmp {output.cm}  # Move uncalibrated CM to final output
            exit 0
        fi

        cat "$estimate_output" >> {log}

        predicted_mem_mb=$(awk '/{params.model_name}/ {{ print $NF; exit }}' "$estimate_output")
        if [[ -z "$predicted_mem_mb" ]]; then
            echo "Could not parse memory estimate; leaving model uncalibrated." >> {log}
            cmstat {input.cm} >> {output.cm_stats} 2>&1 | tee -a {log}
            mv {output.cm}.tmp {output.cm}  # Move uncalibrated CM to final output
            exit 0
        fi
        required_mem_mb=$(awk -v m="$predicted_mem_mb" -v f="{resources.mem_safety_factor}" 'BEGIN {{ printf "%.0f", m * f }}')

        echo "Predicted memory requirement ({threads} cpus, x{resources.mem_safety_factor} safety factor): ${{required_mem_mb}} Mb (budget: {resources.mem_mb} Mb)" >> {log}

        if (( required_mem_mb > {resources.mem_mb} )); then
            echo "Memory requirement too high: ${{required_mem_mb}} Mb > {resources.mem_mb} Mb budget. Leaving uncalibrated." >> {log}
            cmstat {input.cm} >> {output.cm_stats} 2>&1 | tee -a {log}
            mv {output.cm}.tmp {output.cm}  # Move uncalibrated CM to final output
            exit 0
        fi

        predicted_runtime=$(awk '/{params.model_name}/ && $2 ~ /^[0-9]+:[0-9][0-9]:[0-9][0-9]$/ {{ print $2; exit }}' "$estimate_output")

        if [[ -z "$predicted_runtime" ]]; then
            echo "Could not parse runtime forecast; proceeding because memory requirement is within budget." >> {log}
            predicted_seconds=0
        else
            IFS=: read -r h m s <<< "$predicted_runtime"
            predicted_seconds=$((10#$h * 3600 + 10#$m * 60 + 10#$s))
        fi

        if (( predicted_seconds > {params.max_est_runtime} )); then
            echo "Runtime too long: $predicted_runtime leaving uncalibrated." >> {log}
            cmstat {input.cm} >> {output.cm_stats} 2>&1 | tee -a {log}
            mv {output.cm}.tmp {output.cm}  # Move uncalibrated CM to final output
            
        else
            echo "Runtime acceptable: $predicted_runtime" >> {log}

            # Calibrate with Infernal 1.1
            echo "Running cmcalibrate with {threads} threads..." >> {log}
            cmcalibrate --cpu {threads} --tailp {params.tailp} -L {params.total_length_random_seq} {output.cm}.tmp 2>&1 | tee -a {log}
            # Show calibrated stats
            echo "" >> {log}
            echo "Calibrated CM statistics:" >> {log}
            cmstat {output.cm}.tmp >> {output.cm_stats} 2>&1 | tee -a {log}
            echo "Calibration complete for {wildcards.cluster_id}" >> {log}
        fi
        mv {output.cm}.tmp {output.cm}  # Move uncalibrated CM to final output
        """

rule build_final_report:
    input:
        rscape_dirs = lambda wildcards: expand(
            "results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/03_structure/cluster_{cluster_id}_rscape_output",
            target_domain=wildcards.target_domain, region=wildcards.region,
            cluster_id=get_motif_cluster_ids(wildcards)),
        cms = lambda wildcards: expand(
            "results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/04_cms/{cluster_id}_calibrated.cm",
            target_domain=wildcards.target_domain, region=wildcards.region,
            cluster_id=get_covarying_cluster_ids(wildcards)),
        cm_stats = lambda wildcards: expand(
            "results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/04_cms/{cluster_id}_calibrated.cm.stats",
            target_domain=wildcards.target_domain, region=wildcards.region,
            cluster_id=get_covarying_cluster_ids(wildcards))
    output:
        html_report = "results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/final_report.html",
        tsv_report = "results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/candidate_covariations.tsv",
        combined_cms = "results/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/combined_calibrated_cms.cm"
    params:
        script = "workflow/src/build_final_report.py"
    conda:
        "../envs/rna_motif_env.yaml"
    log:
        "logs/03_manual_motif_discovery/{target_domain}_domain_flanking_{region}/final_report.log"
    shell:
        """
        cat {input.cms} > {output.combined_cms} 2> {log}

        if [ -z "{input.rscape_dirs}" ]; then
            echo "No clusters passed filtering (only singletons) - skipping report script." >> {log}
            echo "<!DOCTYPE html><html><head><meta charset=\"utf-8\"><title>{wildcards.target_domain} / {wildcards.region}</title></head><body><h1>{wildcards.target_domain} / {wildcards.region}</h1><p>No candidate clusters found. </p></body></html>" > {output.html_report}
            printf "target_domain\\tregion\\tcluster_id\\tin_structure\\tleft_pos\\tright_pos\\tscore\\tevalue\\tsubstitutions\\tpower\\n" > {output.tsv_report}
        else
            python {params.script} \
                --rscape-dirs {input.rscape_dirs} \
                --cm-stats {input.cm_stats} \
                --target-domain {wildcards.target_domain} --region {wildcards.region} \
                --html-out {output.html_report} \
                --tsv-out {output.tsv_report} \
                2>&1 | tee -a {log}
        fi
        """