"""
RNA structural motif discovery using CMfinder
Based on Narunsky et al., NAR 2024 (doi: 10.1093/nar/gkae248)
"""

# Checkpoint: Initial CMfinder run to identify candidate RNA structures
checkpoint cmfinder_initial:
    input:
        fasta = "results/02_sequence_selection/03_combined/{target_domain}/{target_domain}_domain_flanking_upstream_filtered.fasta"
    output:
        outdir = directory("results/03_motif_discovery/{target_domain}/01_initial_discovery"),
        flag = touch("results/03_motif_discovery/{target_domain}/01_initial_discovery/.complete")
    params:
        basename = "{target_domain}"
    threads: 4
    resources:
        mem_mb=8000,
        runtime=240
    log:
        "logs/cmfinder/{target_domain}_01_initial.log"
    shell:
        """
        mkdir -p {output.outdir}
        
        # Copy input file
        cp {input.fasta} {output.outdir}/{params.basename}.fasta
        
        # Run CMfinder - generates .cm.h* and .motif.h* files
        docker run --rm --platform linux/amd64 \
            -v "$(pwd)/{output.outdir}":/data \
            --entrypoint /bin/bash \
            cmfinder:latest \
            -c "cd /data && /opt/cmfinder-0.4.1.18/bin/cmfinder.pl {params.basename}.fasta" \
            2>&1 | tee {log}
        
        echo "" >> {log}
        echo "CMfinder initial run completed. Generated files:" >> {log}
        ls -lh {output.outdir}/*.cm.h* {output.outdir}/*.motif.h* 2>/dev/null >> {log} || echo "No motif files found" >> {log}
        
        """

# Helper function to get motif IDs from checkpoint metadata
def get_motifs(wildcards):
    import json
    import re
    from pathlib import Path

    motif_dir = Path(checkpoints.cmfinder_initial.get(target_domain=wildcards.target_domain).output.outdir)
    cm_files = list(motif_dir.glob("*.cm.h*"))

    print(cm_files)

    # cm file will be named like: {target_domain}.fasta.cm.h1_1, {target_domain}.fasta.cm.h2_1, etc., ids we want are the h1_1, h2_1, etc. parts

    motif_ids = []
    for cm_file in cm_files:
        match = re.search(r"\.cm\.(h\d+_\d+)$", cm_file.name)
        if match:
            motif_ids.append(match.group(1))

    motif_ids = sorted(set(motif_ids))

    print(f"Identified motif IDs for {wildcards.target_domain}: {motif_ids}")

    return motif_ids

# Rule 3: Calibrate individual covariance models
rule calibrate_cm_model:
    input:
        cm = "results/03_motif_discovery/{target_domain}/01_initial_discovery/{target_domain}.fasta.cm.{motif_id}"
    output:
        cm = "results/03_motif_discovery/{target_domain}/02_models/calibrated/{motif_id}.cm"
    conda:
        "../envs/rna_motif_env.yaml"
    threads: 4
    resources:
        mem_mb=8000,
        runtime=120
    log:
        "logs/cmfinder/{target_domain}_03_calibrate_{motif_id}.log"
    shell: # NEED TO BUILD CM FROM .MOTIF FILES, CMFINDER MAKES TOO OLD FORMAT!!!!
        """
        echo "Calibrating motif {wildcards.motif_id}" > {log}
        
        # convert to modern Infernal format if needed (some older CMfinder versions produce older format)
        cmconvert -o {output.cm}.converted {input.cm}

        # Calibrate in place with  Infernal
        cmcalibrate --cpu {threads} {output.cm}.converted 2>&1 | tee -a {log}
        
        echo "Calibration complete for {wildcards.motif_id}" >> {log}
        """

# Rule 4: Search for homologs with individual models
rule search_for_homologs:
    input:
        cm = "results/03_motif_discovery/{target_domain}/02_models/calibrated/{motif_id}.cm",
        database = "results/02_sequence_selection/03_combined/{target_domain}/{target_domain}_domain_flanking_upstream_filtered.fasta"
    output:
        hits_tbl = "results/03_motif_discovery/{target_domain}/03_homolog_search/hits/{motif_id}_hits.tblout",
        hits_sto = "results/03_motif_discovery/{target_domain}/03_homolog_search/hits/{motif_id}_hits.sto"
    params:
        evalue = lambda wildcards: config.get("cmfinder", {}).get("evalue_threshold", 0.01),
        calibrated_dir = "results/03_motif_discovery/{target_domain}/02_models/calibrated"
    threads: 4
    resources:
        mem_mb=8000,
        runtime=120
    log:
        "logs/cmfinder/{target_domain}_04_cmsearch_{motif_id}.log"
    shell:
        """
        echo "Searching for homologs of motif {wildcards.motif_id}" > {log}
        echo "E-value threshold: {params.evalue}" >> {log}
        
        outdir=$(dirname {output.hits_tbl})
        mkdir -p $outdir
        
        docker run --rm --platform linux/amd64 \
            -v "$(pwd)/$outdir":/output \
            -v "$(pwd)/{params.calibrated_dir}":/models:ro \
            -v "$(pwd)/{input.database}":/data/database.fasta:ro \
            cmfinder:latest \
            cmsearch \
            --cpu {threads} \
            -E {params.evalue} \
            --tblout /output/{wildcards.motif_id}_hits.tblout \
            -A /output/{wildcards.motif_id}_hits.sto \
            /models/{wildcards.motif_id}.cm \
            /data/database.fasta \
            2>&1 | tee -a {log}
        
        # Count hits
        hits=$(grep -v "^#" {output.hits_tbl} | wc -l)
        echo "Found $hits hits for {wildcards.motif_id}" >> {log}
        echo "Found $hits hits for {wildcards.motif_id}"
        """

# Rule 5: Extract homolog sequences
rule extract_homolog_sequences:
    input:
        hits_tbl = "results/03_motif_discovery/{target_domain}/03_homolog_search/hits/{motif_id}_hits.tblout",
        database = "results/02_sequence_selection/03_combined/{target_domain}/{target_domain}_domain_flanking_upstream_filtered.fasta"
    output:
        fasta = "results/03_motif_discovery/{target_domain}/03_homolog_search/expanded_seqs/{motif_id}_expanded.fasta",
        ids = "results/03_motif_discovery/{target_domain}/03_homolog_search/expanded_seqs/{motif_id}_ids.json"
    params:
        script = "workflow/src/extract_motif_sequences.py"
    log:
        "logs/cmfinder/{target_domain}_05_extract_{motif_id}.log"
    conda:
        "../envs/domain_analysis_env.yaml"
    shell:
        """
        python {params.script} \
        --fasta-file {input.database} \
        --motif-hits-tbl {input.hits_tbl} \
        --output-fasta {output.fasta} \
        --output-ids {output.ids} \
        2>&1 | tee {log}
        """

# Rule 6: Refine alignment with expanded sequences
rule refine_alignment:
    input:
        fasta = "results/03_motif_discovery/{target_domain}/03_homolog_search/expanded_seqs/{motif_id}_expanded.fasta"
    output:
        outdir = directory("results/03_motif_discovery/{target_domain}/04_refinement/{motif_id}"),
        flag = touch("results/03_motif_discovery/{target_domain}/04_refinement/{motif_id}/.complete")
    threads: 4
    resources:
        mem_mb=16000,
        runtime=360
    log:
        "logs/cmfinder/{target_domain}_06_refine_{motif_id}.log"
    shell:
        """
        mkdir -p {output.outdir}
        
        # Check if we have sequences to refine
        seq_count=$(grep -c "^>" {input.fasta} || echo 0)
        
        if [ $seq_count -lt 3 ]; then
            echo "Not enough sequences ($seq_count) for refinement of {wildcards.motif_id}" > {log}
            echo "Skipping refinement step" >> {log}
            touch {output.outdir}/insufficient_sequences.txt
            exit 0
        fi
        
        echo "Refining alignment for motif {wildcards.motif_id}" > {log}
        echo "Sequence count: $seq_count" >> {log}
        
        # Copy to output directory
        cp {input.fasta} {output.outdir}/{wildcards.motif_id}.fasta
        
        # Run CMfinder again for refinement
        docker run --rm --platform linux/amd64 \
            -v "$(pwd)/{output.outdir}":/data \
            --entrypoint /bin/bash \
            cmfinder:latest \
            -c "cd /data && /opt/cmfinder-0.4.1.18/bin/cmfinder.pl {wildcards.motif_id}.fasta" \
            2>&1 | tee -a {log}
        
        echo "Refinement complete for {wildcards.motif_id}" >> {log}
        ls -lh {output.outdir}/ >> {log}
        """

# Rule 7: Build final covariance model
rule build_final_model:
    input:
        refined_dir = "results/03_motif_discovery/{target_domain}/04_refinement/{motif_id}",
        flag = "results/03_motif_discovery/{target_domain}/04_refinement/{motif_id}/.complete"
    output:
        cm = "results/03_motif_discovery/{target_domain}/02_models/final/{motif_id}_final.cm"
    log:
        "logs/cmfinder/{target_domain}_07_final_{motif_id}.log"
    shell:
        """
        echo "Building final model for motif {wildcards.motif_id}" > {log}
        
        outdir=$(dirname {output.cm})
        mkdir -p $outdir
        
        # Check if refinement produced a Stockholm file
        sto_file=$(find {input.refined_dir} -name "*.sto" | head -1)
        
        if [ -z "$sto_file" ] || [ ! -f "$sto_file" ]; then
            echo "No Stockholm alignment found, using original model" >> {log}
            # Copy the calibrated model as final model
            cp results/03_motif_discovery/{wildcards.target_domain}/02_models/calibrated/{wildcards.motif_id}.cm {output.cm}
        else
            echo "Building CM from refined alignment: $sto_file" >> {log}
            
            docker run --rm --platform linux/amd64 \
                -v "$(pwd)/{input.refined_dir}":/input:ro \
                -v "$(pwd)/$outdir":/output \
                cmfinder:latest \
                cmbuild /output/{wildcards.motif_id}_final.cm /input/$(basename $sto_file) \
                2>&1 | tee -a {log}
        fi
        
        echo "Final model complete for {wildcards.motif_id}" >> {log}
        """

# Rule 8: Visualize structure
rule visualize_structure:
    input:
        refined_dir = "results/03_motif_discovery/{target_domain}/04_refinement/{motif_id}",
        flag = "results/03_motif_discovery/{target_domain}/04_refinement/{motif_id}/.complete"
    output:
        pdf = "results/03_motif_discovery/{target_domain}/05_visualizations/{motif_id}_structure.pdf"
    log:
        "logs/cmfinder/{target_domain}_08_visualize_{motif_id}.log"
    shell:
        """
        echo "Generating visualization for motif {wildcards.motif_id}" > {log}
        
        outdir=$(dirname {output.pdf})
        mkdir -p $outdir
        
        # Find Stockholm file
        sto_file=$(find {input.refined_dir} -name "*.sto" | head -1)
        
        if [ -z "$sto_file" ] || [ ! -f "$sto_file" ]; then
            echo "No Stockholm alignment found for visualization" >> {log}
            # Create a placeholder
            touch {output.pdf}
        else
            echo "Visualizing: $sto_file" >> {log}
            
            docker run --rm --platform linux/amd64 \
                -v "$(pwd)/{input.refined_dir}":/input:ro \
                -v "$(pwd)/$outdir":/output \
                cmfinder:latest \
                r2r --GSC-weighted-consensus \
                /input/$(basename $sto_file) \
                /output/{wildcards.motif_id}_structure.pdf \
                2>&1 | tee -a {log} || echo "R2R visualization failed" >> {log}
        fi
        """

# Rule 9: Gather all motif analysis outputs (triggers parallel execution)
rule gather_motif_analysis:
    input:
        final_models = lambda wildcards: expand(
            "results/03_motif_discovery/{target_domain}/02_models/final/{motif_id}_final.cm",
            target_domain=wildcards.target_domain,
            motif_id=get_motifs(wildcards)
        ),
        visualizations = lambda wildcards: expand(
            "results/03_motif_discovery/{target_domain}/05_visualizations/{motif_id}_structure.pdf",
            target_domain=wildcards.target_domain,
            motif_id=get_motifs(wildcards)
        ),
        search_results = lambda wildcards: expand(
            "results/03_motif_discovery/{target_domain}/03_homolog_search/hits/{motif_id}_hits.tblout",
            target_domain=wildcards.target_domain,
            motif_id=get_motifs(wildcards)
        )
    output:
        touch("results/03_motif_discovery/{target_domain}/.motif_analysis_complete")
    log:
        "logs/cmfinder/{target_domain}_09_gather.log"
    shell:
        """
        echo "All motif analysis complete for {wildcards.target_domain}" > {log}
        echo "Final models: $(echo {input.final_models} | wc -w)" >> {log}
        echo "Visualizations: $(echo {input.visualizations} | wc -w)" >> {log}
        echo "Search results: $(echo {input.search_results} | wc -w)" >> {log}
        """

# Rule 10: Summarize all results
rule summarize_cmfinder_results:
    input:
        gather_flag = "results/03_motif_discovery/{target_domain}/.motif_analysis_complete"
    output:
        summary = "results/03_motif_discovery/{target_domain}/06_reports/analysis_summary.txt",
        report = "results/03_motif_discovery/{target_domain}/06_reports/report.html"
    params:
        script = "workflow/src/summarize_cmfinder.py",
        cmfinder_dir = "results/03_motif_discovery/{target_domain}"
    log:
        "logs/cmfinder/{target_domain}_10_summary.log"
    conda:
        "../envs/rna_motif_env.yaml"
    shell:
        """
        python {params.script} \
            --target-domain {wildcards.target_domain} \
            --cmfinder-dir {params.cmfinder_dir} \
            --output-summary {output.summary} \
            --output-report {output.report} \
            2>&1 | tee {log}
        """