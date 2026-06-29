"""
RNA structural motif discovery using CMfinder
Based on Narunsky et al., NAR 2024 (doi: 10.1093/nar/gkae248)
"""

# Checkpoint: Initial CMfinder run to identify candidate RNA structures
checkpoint cmfinder_initial:
    input:
        fasta = "results/02_sequence_selection/03_combined/{target_domain}/{target_domain}_{region}_filtered.fasta"
    output:
        outdir = directory("results/03_motif_discovery/{target_domain}_{region}/01_initial_discovery/aligmnments"),
        motif_list = "results/03_motif_discovery/{target_domain}_{region}/01_initial_discovery/motif_list.txt"
    params:
        basename = "{target_domain}_{region}",
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

    # Access checkpoint output
    checkpoint_output = checkpoints.cmfinder_initial.get(**wildcards).output
    motif_list_file = Path(checkpoint_output.motif_list)
    
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
        msa = "results/03_motif_discovery/{target_domain}_{region}/01_initial_discovery/alignments/{initial_motif}",
        flag = "results/03_motif_discovery/{target_domain}_{region}/01_initial_discovery/.complete"
    output:
        cm = "results/03_motif_discovery/{target_domain}_{region}/01_initial_discovery/models/{initial_motif}.cm",
    conda:
        "../envs/rna_motif_env.yaml"
    log:
        "logs/cmfinder/{target_domain}_{region}_02_build_cm_{initial_motif}.log"
    shell:
        """
        echo "Building CM with cmbuild..." >> {log}
        cmbuild -F {output.cm} {input.msa} 2>&1 | tee -a {log}
        echo "CM building complete for {wildcards.initial_motif}" >> {log} 
        """

# Rule 3: Calibrate covariance models (OPTIONAL - can skip for faster results)
rule calibrate_cm_model:
    input:
        cm = "results/03_motif_discovery/{target_domain}_{region}/01_initial_discovery/models/{initial_motif}.cm"
    output:
        cm = "results/03_motif_discovery/{target_domain}_{region}/01_initial_discovery/models/{initial_motif}_calibrated.cm"
    conda:
        "../envs/rna_motif_env.yaml"
    threads: config["cmfinder"]["threads"]["calibration"]
    resources:
        mem_mb= config["cmfinder"]["resources"]["calibration_mb"],
        runtime= config["cmfinder"]["resources"]["calibration_runtime"]
    log:
        "logs/cmfinder/{target_domain}_{region}_03_calibrate_{initial_motif}.log"
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
        cm = "results/03_motif_discovery/{target_domain}_{region}/01_initial_discovery/models/{initial_motif}_calibrated.cm",
        database = "results/02_sequence_selection/03_combined/{target_domain}/{target_domain}_{region}_filtered.fasta"
    output:
        hits_tbl = "results/03_motif_discovery/{target_domain}_{region}/02_homolog_search/hits/{initial_motif}_hits.tblout",
        hits_sto = "results/03_motif_discovery/{target_domain}_{region}/02_homolog_search/hits/{initial_motif}_hits.sto"
    params:
        evalue = config["cmfinder"]["evalue_threshold"]
    conda:
        "../envs/rna_motif_env.yaml"
    threads: config["cmfinder"]["threads"]["homolog_search"]
    resources:
        mem_mb= config["cmfinder"]["resources"]["homolog_search_mb"],
        runtime= config["cmfinder"]["resources"]["homolog_search_runtime"]
    log:
        "logs/cmfinder/{target_domain}_{region}_04_cmsearch_{initial_motif}.log"
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
        hits_tbl = "results/03_motif_discovery/{target_domain}_{region}/02_homolog_search/hits/{initial_motif}_hits.tblout",
        database = "results/02_sequence_selection/03_combined/{target_domain}/{target_domain}_{region}_filtered.fasta"
    output:
        fasta = "results/03_motif_discovery/{target_domain}_{region}/02_homolog_search/expanded_seqs/{initial_motif}_expanded.fasta",
        ids = "results/03_motif_discovery/{target_domain}_{region}/02_homolog_search/expanded_seqs/{initial_motif}_ids.txt"
    conda:
        "../envs/rna_motif_env.yaml"
    params:
        script = "workflow/src/extract_motif_sequences.py"
    log:
        "logs/cmfinder/{target_domain}_{region}_05_extract_{initial_motif}.log"
    shell:
        """
        python {params.script} \
            --fasta-file {input.database} \
            --mmotif-hits-tbl {input.hits_tbl} \
            --output-fasta {output.fasta} \
            --output-ids {output.ids} \
            2>&1 | tee {log}"""
        # """
        # mkdir -p $(dirname {output.fasta})
        
        # echo "Extracting sequences for motif {wildcards.initial_motif}" > {log}
        
        # # Extract sequence IDs from hits (skip header lines)
        # grep -v "^#" {input.hits_tbl} | awk '{{print $1}}' | sort -u > {output.ids}
        
        # count=$(wc -l < {output.ids})
        # echo "Extracted $count unique sequence IDs" >> {log}
        
        # if [ $count -gt 0 ]; then
        #     # Extract sequences using grep
        #     > {output.fasta}
        #     while read seq_id; do
        #         grep -A 1 "^>$seq_id" {input.database} >> {output.fasta}
        #     done < {output.ids}
            
        #     final_count=$(grep -c "^>" {output.fasta})
        #     echo "Extracted $final_count sequences to {output.fasta}" >> {log}
        # else
        #     echo "No hits found, creating empty file" >> {log}
        #     touch {output.fasta}
        # fi
        # """

# Rule 6a: Check if refinement is feasible (sequence count threshold)
# rule check_refinement_feasibility:
#     input:
#         fasta = "results/03_motif_discovery/{target_domain}_{region}/03_homolog_search/expanded_seqs/{initial_motif}_expanded.fasta"
#     output:
#         status = "results/03_motif_discovery/{target_domain}_{region}/04_refinement/{initial_motif}/.status"
#     params:
#         min_refinement_sequences = config["cmfinder"]["min_refinement_sequences"]
#     conda:
#         "../envs/rna_motif_env.yaml"
#     log:
#         "logs/cmfinder/{target_domain}_{region}_06a_check_feasibility_{initial_motif}.log"
#     shell:
#         """
#         seq_count=$(grep -c "^>" {input.fasta})
        
#         echo "Checking refinement feasibility for motif {wildcards.initial_motif}" > {log}
#         echo "Sequence count: $seq_count" >> {log}
#         echo "Minimum required: {params.min_refinement_sequences}" >> {log}
        
#         if [ $seq_count -lt {params.min_refinement_sequences} ]; then
#             echo "SKIP: Not enough sequences for refinement" >> {log}
#             echo "SKIP" > {output.status}
#         else
#             echo "PROCEED: Sufficient sequences for refinement" >> {log}
#             echo "PROCEED" > {output.status}
#         fi
#         """

# Rule 6b: Refine alignment with CMfinder using manifest
checkpoint refine_alignment_cmfinder:
    input:
        fasta = lambda wildcards: expand(
            "results/03_motif_discovery/{target_domain}_{region}/03_homolog_search/expanded_seqs/{initial_motif}_expanded.fasta",
            initial_motif=get_initial_motifs(wildcards),
            target_domain=wildcards.target_domain,
            region=wildcards.region
        )
            #"results/03_motif_discovery/{target_domain}_{region}/03_homolog_search/expanded_seqs/{initial_motif}_expanded.fasta"
    output:
        outdir = directory("results/03_motif_discovery/{target_domain}_{region}/04_refinement"),
        motif_list = "results/03_motif_discovery/{target_domain}_{region}/04_refinement/refined_motif_list.txt"
    params:
        docker_image = config["cmfinder"]["docker_image"],
        docker_platform = config["cmfinder"]["docker_platform"]
    threads: config["cmfinder"]["threads"]["refinement"]
    resources:
        mem_mb=config["cmfinder"]["resources"]["refinement_mb"],
        runtime=config["cmfinder"]["resources"]["refinement_runtime"]
    log:
        "logs/cmfinder/{target_domain}_{region}_06b_refine_cmfinder_{initial_motif}.log"
    shell:
        """
        status=$(cat {input.status})
        echo "Refinement status: $status" > {log}
        
        if [ "$status" = "SKIP" ]; then
            echo "Refinement skipped due to insufficient sequences" >> {log}
            touch {output.outdir}/insufficient_sequences.txt
            touch {output.motif_list}
            touch {output.status_flag}
            exit 0
        fi
        
        echo "Running CMfinder refinement for motif {wildcards.initial_motif}" >> {log}
        
        # Copy to output directory
        cp {input.fasta} {output.outdir}/{wildcards.initial_motif}.fasta
        
        # Run CMfinder refinement with higher stringency and motif list output
        docker run --rm --platform {params.docker_platform} \
            -v "$(pwd)/{output.outdir}":/data \
            --entrypoint /bin/bash \
            {params.docker_image} \
            -c "cd /data && /opt/cmfinder-0.4.1.18/bin/cmfinder04.pl -f 0.8 -cpu {threads} -motifList refined_motif_list.txt {wildcards.initial_motif}.fasta" \
            2>&1 | tee -a {log}
        
        echo "Refinement complete for {wildcards.initial_motif}" >> {log}
        echo "Generated motif list:" >> {log}
        cat {output.outdir}/refined_motif_list.txt >> {log}
        """

# Helper function to get refined motifs
def get_refined_motifs(wildcards):
    from pathlib import Path
    
    checkpoint_output = checkpoints.refine_alignment_cmfinder.get(**wildcards).output
    motif_list_file = Path(checkpoint_output.motif_list)
    
    refined_motif_ids = []
    
    # If refinement was skipped, return empty list
    insufficient_flag = Path(checkpoint_output.outdir) / "insufficient_sequences.txt"
    if insufficient_flag.exists():
        print(f"Refinement was skipped for {wildcards.initial_motif} due to insufficient sequences")
        return []
    
    # Read refined motif list
    if motif_list_file.exists():
        with open(motif_list_file) as f:
            for line in f:
                line = line.strip()
                refined_motif_ids.append(line)
    
    refined_motif_ids = sorted(set(refined_motif_ids))
    print(f"Loaded {len(refined_motif_ids)} refined motifs for {wildcards.initial_motif}: {refined_motif_ids}")
    
    return refined_motif_ids

# Rule 7: Build final covariance model from refined alignment or use calibrated model
rule build_final_model:
    input:
        refined_dir = "results/03_motif_discovery/{target_domain}_{region}/04_refinement/{initial_motif}/{refined_motif}",
        refinement_flag = "results/03_motif_discovery/{target_domain}_{region}/04_refinement/{initial_motif}/.refinement_complete",
        calibrated_cm = "results/03_motif_discovery/{target_domain}_{region}/02_models/calibrated/{initial_motif}.cm"
    output:
        cm = "results/03_motif_discovery/{target_domain}_{region}/04_refinement/final/{motif_id}_final.cm"
    threads: config["cmfinder"]["threads"]["calibration"]
    conda:
        "../envs/rna_motif_env.yaml"
    log:
        "logs/cmfinder/{target_domain}_{region}_07_final_{motif_id}.log"
    shell:
        """
        echo "Building final model for motif {wildcards.motif_id}" > {log}
        
        # Check if refinement was skipped
        if [ -f {input.refined_dir}/insufficient_sequences.txt ]; then
            echo "Refinement was skipped, using calibrated model as final" >> {log}
            cp {input.calibrated_cm} {output.cm}
            echo "Using pre-calibrated model" >> {log}
            cmstat {output.cm} >> {log}
            exit 0
        fi
        
        # Look for refined .motif files (Stockholm alignments) from refined_motif_list.txt
        motif_list="{input.refined_dir}/refined_motif_list.txt"
        
        if [ ! -f "$motif_list" ] || [ ! -s "$motif_list" ]; then
            echo "No refined motif list found, using calibrated model as final" >> {log}
            cp {input.calibrated_cm} {output.cm}
            cmstat {output.cm} >> {log}
            exit 0
        fi
        
        # Extract first refined motif from list
        refined_motif=$(head -1 "$motif_list" | tr -d '\n' | sed 's/[[:space:]]*$//')
        refined_motif_path="{input.refined_dir}/$refined_motif"
        
        if [ -z "$refined_motif" ] || [ ! -f "$refined_motif_path" ]; then
            echo "Refined motif file not found: $refined_motif_path" >> {log}
            echo "Using calibrated model as fallback" >> {log}
            cp {input.calibrated_cm} {output.cm}
            cmstat {output.cm} >> {log}
            exit 0
        fi
        
        echo "Building final CM from refined alignment: $refined_motif" >> {log}
        
        # Build and calibrate final CM using Infernal
        cmbuild -F {output.cm} "$refined_motif_path" 2>&1 | tee -a {log}
        
        echo "Calibrating final model..." >> {log}
        cmcalibrate --cpu {threads} {output.cm} 2>&1 | tee -a {log}
        
        echo "Final model complete for {wildcards.motif_id}" >> {log}
        cmstat {output.cm} >> {log}
        """

# Rule 8: Visualize structure with R2R (if available)
# Rule 8: Visualize structure with multiple fallback options
rule visualize_structure:
    input:
        refined_dir = "results/03_motif_discovery/{target_domain}_{region}/04_refinement/{motif_id}",
        #original_alignment = "results/03_motif_discovery/{target_domain}_{region}/02_models/built/{motif_id}.sto",
        refinement_flag = "results/03_motif_discovery/{target_domain}_{region}/04_refinement/{motif_id}/.refinement_complete"
    output:
        pdf = "results/03_motif_discovery/{target_domain}_{region}/05_visualizations/{motif_id}_structure.pdf",
        svg = "results/03_motif_discovery/{target_domain}_{region}/05_visualizations/{motif_id}_structure.svg",
        txt = "results/03_motif_discovery/{target_domain}_{region}/05_visualizations/{motif_id}_structure.txt"
    conda:
        "../envs/rna_motif_env.yaml"
    params:
        script = "workflow/src/visualize_structure.py"
    log:
        "logs/cmfinder/{target_domain}_{region}_08_visualize_{motif_id}.log"
    shell:
        """
        echo "Generating visualization for motif {wildcards.motif_id}" > {log}
        
        # Determine which alignment to use
        if [ -f {input.refined_dir}/insufficient_sequences.txt ]; then
            echo "Refinement was skipped, using original alignment" >> {log}
            sto_file="{input.original_alignment}"
        else
            # Find refined Stockholm alignment
            sto_file=$(find {input.refined_dir} -name "*.motif.h*" | head -1)
            if [ -z "$sto_file" ] || [ ! -f "$sto_file" ]; then
                echo "No refined alignment, using original" >> {log}
                sto_file="{input.original_alignment}"
            fi
        fi
        
        echo "Using alignment: $sto_file" >> {log}
        
        # Method 1: Try R2R (often crashes)
        echo "Attempting R2R visualization..." >> {log}
        if command -v r2r &> /dev/null; then
            if r2r --GSC-weighted-consensus "$sto_file" {output.pdf} 2>> {log}; then
                echo "R2R visualization successful" >> {log}
                # Also try SVG
                r2r --GSC-weighted-consensus "$sto_file" {output.svg} 2>> {log} || echo "R2R SVG failed" >> {log}
            else
                echo "R2R crashed or failed" >> {log}
            fi
        else
            echo "R2R not available" >> {log}
        fi
        
        # Method 2: Use Infernal's esl-alipid for text visualization
        echo "" >> {log}
        echo "Creating text-based structure visualization..." >> {log}
        if command -v esl-alipid &> /dev/null; then
            esl-alipid "$sto_file" > {output.txt} 2>> {log}
            echo "Text visualization created: {output.txt}" >> {log}
        else
            echo "esl-alipid not available, extracting structure from Stockholm file" >> {log}
            # Extract structure annotation from Stockholm file
            grep -E "^#=GC SS_cons|^#=GR .* SS" "$sto_file" > {output.txt} 2>> {log} || \
                echo "Could not extract structure" > {output.txt}
        fi
        
        # Method 3: Create a simple PDF from text if R2R failed
        if [ ! -s {output.pdf} ]; then
            echo "Creating fallback PDF visualization..." >> {log}

            python {params.script} \
                --stockholm "$sto_file" \
                --output {output.pdf} \
                --motif-id {wildcards.motif_id} \
                2>> {log}

            echo "Fallback PDF visualization created" >> {log}
        fi
        
        # Ensure all output files exist
        touch {output.pdf} {output.svg} {output.txt}
        
        echo "Visualization complete" >> {log}
        """



rule annotate_motif_families_comprehensive:
    input:
        expanded_seqs = "results/03_motif_discovery/{target_domain}_{region}/03_homolog_search/expanded_seqs/{motif_id}_expanded.fasta",
        original_alignment = "results/03_motif_discovery/{target_domain}_{region}/02_models/built/{motif_id}.sto",
        rfam_db = "resources/rfam/Rfam.cm"
    output:
        rfam_hits = "results/03_motif_discovery/{target_domain}_{region}/06_annotation/{motif_id}_rfam_hits.tblout",
        rfam_summary = "results/03_motif_discovery/{target_domain}_{region}/06_annotation/{motif_id}_rfam_summary.txt",
        novel_flag = "results/03_motif_discovery/{target_domain}_{region}/06_annotation/{motif_id}_novelty.txt"
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
        "logs/cmfinder/{target_domain}_{region}_08_annotate_{motif_id}.log"
    shell:
        """
        echo "Comprehensive Rfam annotation for motif {wildcards.motif_id}" > {log}
        
        # Index Rfam if needed
        if [ ! -f {input.rfam_db}.i1m ]; then
            echo "Indexing Rfam database..." >> {log}
            cmpress {input.rfam_db} 2>&1 | tee -a {log}
        fi
        
        # Decide which sequences to search
        if [ -s {input.expanded_seqs} ]; then
            search_file="{input.expanded_seqs}"
            echo "Using expanded sequences (with homologs)" >> {log}
        else
            # Extract sequences from original alignment
            search_file=$(mktemp --suffix=.fasta)
            echo "Extracting sequences from original alignment" >> {log}
            grep -v "^#" {input.original_alignment} | grep -v "^//" | \
                awk '{{print ">"$1"\\n"$2}}' > $search_file
        fi
        
        seq_count=$(grep -c "^>" $search_file 2>/dev/null || echo 0)
        echo "Searching $seq_count sequences against Rfam" >> {log}
        
        # Search against Rfam
        cmscan \
            --cpu {threads} \
            -E {params.evalue} \
            --tblout {output.rfam_hits} \
            {input.rfam_db} \
            $search_file \
            2>&1 | tee -a {log}
        
        # Use Python script to parse and summarize results
        python {params.script} \
            --rfam-hits {output.rfam_hits} \
            --motif-id {wildcards.motif_id} \
            --seq-count $seq_count \
            --output-summary {output.rfam_summary} \
            --output-flag {output.novel_flag} \
            2>&1 | tee -a {log}
        
        # Clean up temp file if created
        if [ ! -s {input.expanded_seqs} ]; then
            rm -f $search_file
        fi
        """

# Rule 9: Gather all motif analysis outputs
rule gather_motif_analysis:
    input:
        final_models = lambda wildcards: expand(
            "results/03_motif_discovery/{target_domain}_{region}/02_models/final/{motif_id}_final.cm",
            target_domain=wildcards.target_domain,
            region=wildcards.region,
            motif_id=get_motifs(wildcards)
        ),
        visualizations = lambda wildcards: expand(
            "results/03_motif_discovery/{target_domain}_{region}/05_visualizations/{motif_id}_structure.pdf",
            target_domain=wildcards.target_domain,
            region=wildcards.region,
            motif_id=get_motifs(wildcards)
        ),
        search_results = lambda wildcards: expand(
            "results/03_motif_discovery/{target_domain}_{region}/03_homolog_search/hits/{motif_id}_hits.tblout",
            target_domain=wildcards.target_domain,
            region=wildcards.region,
            motif_id=get_motifs(wildcards)
        ),
        annotations = lambda wildcards: expand(
            "results/03_motif_discovery/{target_domain}_{region}/06_annotation/{motif_id}_novelty.txt",
            target_domain=wildcards.target_domain,
            region=wildcards.region,
            motif_id=get_motifs(wildcards)
        )
    output:
        touch("results/03_motif_discovery/{target_domain}_{region}/.motif_analysis_complete")
    log:
        "logs/cmfinder/{target_domain}_{region}_09_gather.log"
    shell:
        """
        echo "All motif analysis complete for {wildcards.target_domain}_{wildcards.region}" > {log}
        echo "Final models: $(echo {input.final_models} | wc -w)" >> {log}
        echo "Visualizations: $(echo {input.visualizations} | wc -w)" >> {log}
        echo "Search results: $(echo {input.search_results} | wc -w)" >> {log}
        """

# Rule 10: Summarize all results
rule summarize_cmfinder_results:
    input:
        gather_flag = "results/03_motif_discovery/{target_domain}_{region}/.motif_analysis_complete"
    output:
        summary = "results/03_motif_discovery/{target_domain}_{region}/07_reports/analysis_summary.txt",
        report = "results/03_motif_discovery/{target_domain}_{region}/07_reports/report.html"
    params:
        script = "workflow/src/summarize_cmfinder.py",
        cmfinder_dir = "results/03_motif_discovery/{target_domain}_{region}"
    conda:
        "../envs/rna_motif_env.yaml"
    log:
        "logs/cmfinder/{target_domain}_{region}_10_summary.log"
    shell:
        """
        python {params.script} \
            --target-domain {wildcards.target_domain} \
            --cmfinder-dir {params.cmfinder_dir} \
            --output-summary {output.summary} \
            --output-report {output.report} \
            2>&1 | tee {log}
        """