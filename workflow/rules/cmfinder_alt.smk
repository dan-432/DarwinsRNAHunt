"""
RNA structural motif discovery using CMfinder
Based on Narunsky et al., NAR 2024 (doi: 10.1093/nar/gkae248)
"""

# Checkpoint: Initial CMfinder run to identify candidate RNA structures
checkpoint cmfinder_initial:
    input:
        fasta = "results/02_sequence_selection/{target_domain}/03_combined/{target_domain}_domain_flanking_{region}_filtered.fasta"
    output:
        outdir = directory("results/03_motif_discovery/{target_domain}_domain_flanking_{region}/01_initial_discovery/alignments"),
        motif_list = "results/03_motif_discovery/{target_domain}_domain_flanking_{region}/01_initial_discovery/alignments/motif_list.txt"
    params:
        basename = "{region}",
        docker_image = config["cmfinder"]["docker_image"],
        docker_platform = config["cmfinder"]["docker_platform"],
        expected_initial_frequency = 0.1
    threads: config["cmfinder"]["threads"]["initial_discovery"]
    resources:
        mem_mb= config["cmfinder"]["resources"]["initial_discovery_mb"],
        runtime= config["cmfinder"]["resources"]["initial_discovery_runtime"]
    log:
        "logs/cmfinder/{target_domain}_{region}_01_initial.log"
    shell:
        """
        mkdir -p {output.outdir}
        # Copy input file
        cp {input.fasta} {output.outdir}/{params.basename}.fasta
        
        # Run CMfinder - generates .motif.h* (Stockholm)
        docker run --rm --platform {params.docker_platform} \
            -v "$(pwd)/{output.outdir}":/data \
            --entrypoint /bin/bash \
            {params.docker_image} \
            -c "cd /data && /opt/cmfinder-0.4.1.18/bin/cmfinder04.pl -f {params.expected_initial_frequency} -cpu {threads} -combine -motifList motif_list.txt {params.basename}.fasta" \
            2>&1 | tee {log}
        
        echo "" >> {log}
        echo "CMfinder initial run completed. Generated motif list:" >> {log}
        cat {output.outdir}/motif_list.txt >> {log}
        """

# Helper function to get motif IDs from checkpoint
def get_initial_motifs(wildcards):
    from pathlib import Path

    print(f"Retrieving initial motifs for {wildcards.target_domain} in region {wildcards.region}")

    # Access checkpoint output
    checkpoint_output = checkpoints.cmfinder_initial.get(**wildcards).output
    motif_list_file = Path(checkpoint_output.motif_list)

    print(f"Motif list file path: {motif_list_file}")
    
    # Read motif list from file
    motif_ids = []
    if motif_list_file.exists():
        with open(motif_list_file) as f:
            for line in f:
                line = line.strip()
                motif_ids.append(line)
    
    motif_ids = sorted(set(motif_ids))
    
    print(f"Loaded {len(motif_ids)} motifs for {wildcards.target_domain}: {motif_ids}")
    
    return motif_ids

# Rule 2: Build modern Infernal CM from Stockholm alignment (.motif file)
rule build_cm_from_motif:
    input:
        msa = "results/03_motif_discovery/{target_domain}_domain_flanking_{region}/01_initial_discovery/alignments/{initial_motif}"
    output:
        cm = "results/03_motif_discovery/{target_domain}_domain_flanking_{region}/01_initial_discovery/models/{initial_motif}_uncalibrated.cm",
    conda:
        "../envs/rna_motif_env.yaml"
    log:
        "logs/cmfinder/{target_domain}_domain_flanking_{region}_02_build_cm_{initial_motif}.log"
    shell:
        """
        echo "Building CM with cmbuild..." >> {log}
        cmbuild -F {output.cm} {input.msa} 2>&1 | tee -a {log}
        echo "CM building complete for {wildcards.initial_motif}" >> {log} 
        """

# Rule 3: Calibrate covariance models (OPTIONAL - can skip for faster results)
rule calibrate_cm_model:
    input:
        cm = "results/03_motif_discovery/{target_domain}_domain_flanking_{region}/01_initial_discovery/models/{initial_motif}_uncalibrated.cm"
    output:
        cm = "results/03_motif_discovery/{target_domain}_domain_flanking_{region}/01_initial_discovery/models/{initial_motif}_calibrated.cm"
    conda:
        "../envs/rna_motif_env.yaml"
    threads: config["cmfinder"]["threads"]["calibration"]
    resources:
        mem_mb= config["cmfinder"]["resources"]["calibration_mb"],
        runtime= config["cmfinder"]["resources"]["calibration_runtime"]
    log:
        "logs/cmfinder/{target_domain}_domain_flanking_{region}_03_calibrate_{initial_motif}.log"
    shell:
        """
        echo "Calibrating motif {wildcards.initial_motif}" > {log}
        
        # Copy CM to calibration directory
        cp {input.cm} {output.cm}
        
        # Calibrate with Infernal 1.1
        echo "Running cmcalibrate with {threads} threads..." >> {log}
        cmcalibrate --cpu {threads} {output.cm} 2>&1 | tee -a {log}
        
        # Show calibrated stats
        echo "" >> {log}
        echo "Calibrated CM statistics:" >> {log}
        cmstat {output.cm} >> {log}
        
        echo "Calibration complete for {wildcards.initial_motif}" >> {log}
        """

# Rule 4: Search for homologs with calibrated models
rule search_for_homologs:
    input:
        cm = "results/03_motif_discovery/{target_domain}_domain_flanking_{region}/01_initial_discovery/models/{initial_motif}_calibrated.cm",
        database = "results/02_sequence_selection/{target_domain}/03_combined/{target_domain}_domain_flanking_{region}_filtered.fasta"
    output:
        hits_tbl = "results/03_motif_discovery/{target_domain}_domain_flanking_{region}/02_homolog_search/hits/{initial_motif}_hits.tblout",
        hits_sto = "results/03_motif_discovery/{target_domain}_domain_flanking_{region}/02_homolog_search/hits/{initial_motif}_hits.sto"
    params:
        evalue = config["cmfinder"]["evalue_threshold"]
    conda:
        "../envs/rna_motif_env.yaml"
    threads: config["cmfinder"]["threads"]["homolog_search"]
    resources:
        mem_mb= config["cmfinder"]["resources"]["homolog_search_mb"],
        runtime= config["cmfinder"]["resources"]["homolog_search_runtime"]
    log:
        "logs/cmfinder/{target_domain}_domain_flanking_{region}_04_cmsearch_{initial_motif}.log"
    shell:
        """
        
        echo "Searching for homologs of motif {wildcards.initial_motif}" > {log}
        echo "E-value threshold: {params.evalue}" >> {log}
        echo "Using {threads} threads" >> {log}
        
        cmsearch \
            --cpu {threads} \
            -E {params.evalue} \
            --tblout {output.hits_tbl} \
            -A {output.hits_sto} \
            {input.cm} \
            {input.database} \
            2>&1 | tee -a {log}
        
        hits=$(grep -v "^#" {output.hits_tbl} | wc -l)
        echo "" >> {log}
        echo "Found $hits hits for {wildcards.initial_motif}" >> {log}
        echo "Found $hits hits for {wildcards.initial_motif}"
        """

# Rule 5: Extract homolog sequences
rule extract_homolog_sequences:
    input:
        hits_tbl = "results/03_motif_discovery/{target_domain}_domain_flanking_{region}/02_homolog_search/hits/{initial_motif}_hits.tblout",
        database = "results/02_sequence_selection/{target_domain}/03_combined/{target_domain}_domain_flanking_{region}_filtered.fasta"
    output:
        fasta = "results/03_motif_discovery/{target_domain}_domain_flanking_{region}/02_homolog_search/expanded_seqs/{initial_motif}_expanded.fasta",
        ids = "results/03_motif_discovery/{target_domain}_domain_flanking_{region}/02_homolog_search/expanded_seqs/{initial_motif}_ids.txt"
    conda:
        "../envs/domain_analysis_env.yaml"
    params:
        script = "workflow/src/extract_motif_sequences.py"
    log:
        "logs/cmfinder/{target_domain}_domain_flanking_{region}_05_extract_{initial_motif}.log"
    shell:
        """
        python {params.script} \
            --fasta-file {input.database} \
            --motif-hits-tbl {input.hits_tbl} \
            --output-fasta {output.fasta} \
            --output-ids {output.ids} \
            2>&1 | tee {log}"""


checkpoint refine_alignment_cmfinder:
    input:
        fasta = "results/03_motif_discovery/{target_domain}_domain_flanking_{region}/02_homolog_search/expanded_seqs/{initial_motif}_expanded.fasta"
    output:
        outdir = directory("results/03_motif_discovery/{target_domain}_domain_flanking_{region}/03_refinement/{initial_motif}"),
        motif_list = "results/03_motif_discovery/{target_domain}_domain_flanking_{region}/03_refinement/{initial_motif}/refined_motif_list.txt"
    params:
        docker_image = config["cmfinder"]["docker_image"],
        docker_platform = config["cmfinder"]["docker_platform"],
        basename = "{initial_motif}_refined",
        min_refinement_sequences = config["cmfinder"]["min_refinement_sequences"]
    threads: config["cmfinder"]["threads"]["refinement"]
    resources:
        mem_mb=config["cmfinder"]["resources"]["refinement_mb"],
        runtime=config["cmfinder"]["resources"]["refinement_runtime"]
    log:
        "logs/cmfinder/{target_domain}_domain_flanking_{region}_06b_refine_{initial_motif}.log"
    shell:
        """
        mkdir -p {output.outdir}
        seq_count=$(grep -c "^>" {input.fasta} || echo 0)

        echo "Refinement for {wildcards.initial_motif}: $seq_count sequences (min {params.min_refinement_sequences})" > {log}

        if [ "$seq_count" -lt "{params.min_refinement_sequences}" ]; then
            echo "Insufficient sequences for refinement, skipping CMfinder" >> {log}
            touch {output.outdir}/insufficient_sequences.txt
            touch {output.motif_list}
        else
            cp {input.fasta} {output.outdir}/{params.basename}.fasta

            docker run --rm --platform {params.docker_platform} \
                -v "$(pwd)/{output.outdir}":/data \
                --entrypoint /bin/bash \
                {params.docker_image} \
                -c "cd /data && /opt/cmfinder-0.4.1.18/bin/cmfinder04.pl -f 0.8 -cpu {threads} -motifList refined_motif_list.txt {params.basename}.fasta" \
                2>&1 | tee -a {log}

            echo "Generated refined motif list:" >> {log}
            cat {output.motif_list} >> {log}
        fi
        """


# Helper: refined motif IDs for one initial_motif (unchanged from original,
# just noting it now matches 6b's nested output structure)
def get_refined_motifs(wildcards):
    from pathlib import Path

    checkpoint_output = checkpoints.refine_alignment_cmfinder.get(**wildcards).output
    outdir = Path(checkpoint_output.outdir)
    motif_list_file = Path(checkpoint_output.motif_list)

    if (outdir / "insufficient_sequences.txt").exists():
        print(f"Refinement skipped for {wildcards.initial_motif} (insufficient sequences)")
        return []

    refined_motif_ids = []
    if motif_list_file.exists():
        with open(motif_list_file) as f:
            refined_motif_ids = sorted(set(line.strip() for line in f if line.strip()))

    print(f"Loaded {len(refined_motif_ids)} refined motifs for {wildcards.initial_motif}: {refined_motif_ids}")
    return refined_motif_ids


# Helper: final CM target paths across ALL initial motifs for a
# (target_domain, region) pair, covering both the "refined" and the
# "skipped -> promoted calibrated model" cases. Used by gather/summary.
def get_final_targets(wildcards):
    from pathlib import Path

    base = f"results/03_motif_discovery/{wildcards.target_domain}_domain_flanking_{wildcards.region}"
    targets = []

    for initial_motif in get_initial_motifs(wildcards):
        print(f"Checking refinement status for initial motif {initial_motif}")
        sub_wc = dict(wildcards)
        sub_wc["initial_motif"] = initial_motif
        checkpoint_output = checkpoints.refine_alignment_cmfinder.get(**sub_wc).output
        outdir = Path(checkpoint_output.outdir)

        if (outdir / "insufficient_sequences.txt").exists():
            targets.append(f"{base}/03_refinement/final/{initial_motif}_final.cm")
        else:
            motif_list_file = Path(checkpoint_output.motif_list)
            if motif_list_file.exists():
                refined_motifs = sorted(set(l.strip() for l in open(motif_list_file) if l.strip()))
                for refined_motif in refined_motifs:
                    # NOTE: refined_motif filenames already contain initial_motif
                    # as a prefix (CMfinder basename = "{initial_motif}_refined"),
                    # so we nest by initial_motif directory rather than
                    # concatenating it into the filename again.
                    targets.append(f"{base}/03_refinement/final/{initial_motif}/{refined_motif}_final.cm")

    return targets


# ---------------------------------------------------------------------
# Rule 7a: Build final CM from a successfully refined alignment
# ---------------------------------------------------------------------
rule build_final_model_refined:
    input:
        refined_msa = "results/03_motif_discovery/{target_domain}_domain_flanking_{region}/03_refinement/{initial_motif}/{refined_motif}",
    output:
        cm = "results/03_motif_discovery/{target_domain}_domain_flanking_{region}/03_refinement/final/{initial_motif}/{refined_motif}_final.cm"
    threads: config["cmfinder"]["threads"]["calibration"]
    conda:
        "../envs/rna_motif_env.yaml"
    log:
        "logs/cmfinder/{target_domain}_domain_flanking_{region}_07a_final_{initial_motif}_{refined_motif}.log"
    shell:
        """
        echo "Building final model from refined alignment {wildcards.refined_motif}" > {log}

        cmbuild -F {output.cm} {input.refined_msa} 2>&1 | tee -a {log}

        echo "Calibrating final model..." >> {log}
        cmcalibrate --cpu {threads} {output.cm} 2>&1 | tee -a {log}

        echo "Final model complete for {wildcards.initial_motif}_{wildcards.refined_motif}" >> {log}
        cmstat {output.cm} >> {log}
        """


# ---------------------------------------------------------------------
# Rule 7b: Refinement was skipped (too few sequences) -> promote the
# already-calibrated initial CM to be the "final" model for this motif.
# ---------------------------------------------------------------------
rule build_final_model_skipped:
    input:
        calibrated_cm = "results/03_motif_discovery/{target_domain}_domain_flanking_{region}/01_initial_discovery/models/{initial_motif}_calibrated.cm",
    output:
        cm = "results/03_motif_discovery/{target_domain}_domain_flanking_{region}/03_refinement/final/{initial_motif}_final.cm"
    conda:
        "../envs/rna_motif_env.yaml"
    log:
        "logs/cmfinder/{target_domain}_domain_flanking_{region}_07b_final_skipped_{initial_motif}.log"
    shell:
        """
        echo "Refinement was skipped for {wildcards.initial_motif}, using calibrated model as final" > {log}
        cp {input.calibrated_cm} {output.cm}
        cmstat {output.cm} >> {log}
        """


# ---------------------------------------------------------------------
# Rule 8: Visualize structure with R2R (falls back to text/basic PDF)
# Matches whichever final-CM pattern produced {motif_id}: split on the
# same refined/skipped distinction so the right source alignment is used.
# ---------------------------------------------------------------------
rule visualize_structure_refined:
    input:
        refined_msa = "results/03_motif_discovery/{target_domain}_domain_flanking_{region}/03_refinement/{initial_motif}/{refined_motif}"
    output:
        pdf = "results/03_motif_discovery/{target_domain}_domain_flanking_{region}/04_visualizations/{initial_motif}/{refined_motif}_structure.pdf"
    conda:
        "../envs/rna_motif_env.yaml"
    params:
        script = "workflow/src/visualise_motifs.py",
        motif_id = "{initial_motif}_{refined_motif}"
    log:
        "logs/cmfinder/{target_domain}_domain_flanking_{region}_08_visualize_{initial_motif}_{refined_motif}.log"
    shell:
        """
        echo "Generating visualization for {params.motif_id}" > {log}

        python {params.script} \
            --stockholm {input.refined_msa} \
            --output {output.pdf} \
            --motif-id {params.motif_id} \
            2>&1 | tee -a {log}

        echo "Visualization complete" >> {log}
        """


rule visualize_structure_skipped:
    input:
        original_alignment = "results/03_motif_discovery/{target_domain}_domain_flanking_{region}/01_initial_discovery/alignments/{initial_motif}"
    output:
        pdf = "results/03_motif_discovery/{target_domain}_domain_flanking_{region}/04_visualizations/{initial_motif}_structure.pdf"
    conda:
        "../envs/rna_motif_env.yaml"
    params:
        script = "workflow/src/visualise_motifs.py"
    log:
        "logs/cmfinder/{target_domain}_domain_flanking_{region}_08_visualize_skipped_{initial_motif}.log"
    shell:
        """
        echo "Generating visualization for {wildcards.initial_motif} (refinement was skipped)" > {log}

        python {params.script} \
            --stockholm {input.original_alignment} \
            --output {output.pdf} \
            --motif-id {wildcards.initial_motif} \
            2>&1 | tee -a {log}

        echo "Visualization complete" >> {log}
        """


# ---------------------------------------------------------------------
# Rule (annotation): comprehensive Rfam annotation, split the same way
# ---------------------------------------------------------------------
rule annotate_motif_refined:
    input:
        expanded_seqs = "results/03_motif_discovery/{target_domain}_domain_flanking_{region}/02_homolog_search/expanded_seqs/{initial_motif}_expanded.fasta",
        rfam_db = "resources/rfam/Rfam.cm"
    output:
        rfam_hits = "results/03_motif_discovery/{target_domain}_domain_flanking_{region}/05_annotation/{initial_motif}/{refined_motif}_rfam_hits.tblout",
        rfam_summary = "results/03_motif_discovery/{target_domain}_domain_flanking_{region}/05_annotation/{initial_motif}/{refined_motif}_rfam_summary.txt",
        novel_flag = "results/03_motif_discovery/{target_domain}_domain_flanking_{region}/05_annotation/{initial_motif}/{refined_motif}_novelty.txt"
    params:
        evalue = 0.01,
        script = "workflow/src/summarize_rfam_hits.py",
        motif_id = "{initial_motif}_{refined_motif}"
    conda:
        "../envs/rna_motif_env.yaml"
    threads: 8
    resources:
        mem_mb=16000,
        runtime=120
    log:
        "logs/cmfinder/{target_domain}_domain_flanking_{region}_annotate_{initial_motif}_{refined_motif}.log"
    shell:
        """
        echo "Rfam annotation for {params.motif_id}" > {log}

        if [ ! -f {input.rfam_db}.i1m ]; then
            cmpress {input.rfam_db} 2>&1 | tee -a {log}
        fi

        seq_count=$(grep -c "^>" {input.expanded_seqs} 2>/dev/null || echo 0)
        echo "Searching $seq_count sequences against Rfam" >> {log}

        cmscan \
            --cpu {threads} \
            -E {params.evalue} \
            --tblout {output.rfam_hits} \
            {input.rfam_db} \
            {input.expanded_seqs} \
            2>&1 | tee -a {log}

        python {params.script} \
            --rfam-hits {output.rfam_hits} \
            --motif-id {params.motif_id} \
            --seq-count $seq_count \
            --output-summary {output.rfam_summary} \
            --output-flag {output.novel_flag} \
            2>&1 | tee -a {log}
        """


rule annotate_motif_skipped:
    input:
        original_alignment = "results/03_motif_discovery/{target_domain}_domain_flanking_{region}/01_initial_discovery/alignments/{initial_motif}",
        rfam_db = "resources/rfam/Rfam.cm"
    output:
        rfam_hits = "results/03_motif_discovery/{target_domain}_domain_flanking_{region}/05_annotation/{initial_motif}_rfam_hits.tblout",
        rfam_summary = "results/03_motif_discovery/{target_domain}_domain_flanking_{region}/05_annotation/{initial_motif}_rfam_summary.txt",
        novel_flag = "results/03_motif_discovery/{target_domain}_domain_flanking_{region}/05_annotation/{initial_motif}_novelty.txt"
    params:
        evalue = 0.01,
        script = "workflow/src/summarize_rfam_hits.py"
    conda:
        "../envs/rna_motif_env.yaml"
    threads: 8
    resources:
        mem_mb=16000,
        runtime=120
    log:
        "logs/cmfinder/{target_domain}_domain_flanking_{region}_annotate_skipped_{initial_motif}.log"
    shell:
        """
        echo "Rfam annotation for {wildcards.initial_motif} (no refinement, using original alignment)" > {log}

        if [ ! -f {input.rfam_db}.i1m ]; then
            cmpress {input.rfam_db} 2>&1 | tee -a {log}
        fi

        search_file=$(mktemp --suffix=.fasta)
        grep -v "^#" {input.original_alignment} | grep -v "^//" | \
            awk '{{print ">"$1"\\n"$2}}' > $search_file
        seq_count=$(grep -c "^>" $search_file 2>/dev/null || echo 0)
        echo "Searching $seq_count sequences against Rfam" >> {log}

        cmscan \
            --cpu {threads} \
            -E {params.evalue} \
            --tblout {output.rfam_hits} \
            {input.rfam_db} \
            $search_file \
            2>&1 | tee -a {log}

        python {params.script} \
            --rfam-hits {output.rfam_hits} \
            --motif-id {wildcards.initial_motif} \
            --seq-count $seq_count \
            --output-summary {output.rfam_summary} \
            --output-flag {output.novel_flag} \
            2>&1 | tee -a {log}

        rm -f $search_file
        """


# ---------------------------------------------------------------------
# Rule 9: Gather all motif analysis outputs
# Uses get_final_targets() to derive matching visualization/annotation
# paths for both refined and skipped motifs, so nothing gets dropped.
# ---------------------------------------------------------------------
def get_visualization_and_annotation_targets(wildcards):
    from pathlib import Path

    finals = get_final_targets(wildcards)
    viz, annot = [], []
    for f in finals:
        p = Path(f)
        base = f.rsplit("/03_refinement/", 1)[0]
        rel_parts = p.relative_to(Path(base) / "03_refinement" / "final").parts

        if len(rel_parts) == 2:
            # Refined case: final/{initial_motif}/{refined_motif}_final.cm
            initial_motif, fname = rel_parts
            refined_motif = fname[: -len("_final.cm")]
            viz.append(f"{base}/04_visualizations/{initial_motif}/{refined_motif}_structure.pdf")
            annot.append(f"{base}/05_annotation/{initial_motif}/{refined_motif}_novelty.txt")
        else:
            # Skipped case: final/{initial_motif}_final.cm
            motif_id = rel_parts[0][: -len("_final.cm")]
            viz.append(f"{base}/04_visualizations/{motif_id}_structure.pdf")
            annot.append(f"{base}/05_annotation/{motif_id}_novelty.txt")

    return viz + annot


rule gather_motif_analysis:
    input:
        final_models = get_final_targets,
        viz_and_annotations = get_visualization_and_annotation_targets
    output:
        touch("results/03_motif_discovery/{target_domain}_domain_flanking_{region}/.motif_analysis_complete")
    log:
        "logs/cmfinder/{target_domain}_domain_flanking_{region}_09_gather.log"
    shell:
        """
        echo "All motif analysis complete for {wildcards.target_domain}_{wildcards.region}" > {log}
        echo "Final models: {input.final_models}" >> {log}
        echo "Visualizations + annotations: {input.viz_and_annotations}" >> {log}
        """


# ---------------------------------------------------------------------
# Rule 10: Summarize all results (unchanged from original)
# ---------------------------------------------------------------------
rule summarize_cmfinder_results:
    input:
        gather_flag = "results/03_motif_discovery/{target_domain}_domain_flanking_{region}/.motif_analysis_complete"
    output:
        summary = "results/03_motif_discovery/{target_domain}_domain_flanking_{region}/06_reports/analysis_summary.txt",
        report = "results/03_motif_discovery/{target_domain}_domain_flanking_{region}/06_reports/report.html"
    params:
        script = "workflow/src/summarize_cmfinder.py",
        cmfinder_dir = "results/03_motif_discovery/{target_domain}_domain_flanking_{region}"
    conda:
        "../envs/rna_motif_env.yaml"
    log:
        "logs/cmfinder/{target_domain}_domain_flanking_{region}_10_summary.log"
    shell:
        """
        python {params.script} \
            --target-domain {wildcards.target_domain} \
            --cmfinder-dir {params.cmfinder_dir} \
            --output-summary {output.summary} \
            --output-report {output.report} \
            2>&1 | tee {log}
        """
