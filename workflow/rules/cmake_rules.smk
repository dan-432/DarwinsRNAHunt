# Complete Snakemake workflow for CMfinder analysis
# Based on Narunsky et al., NAR 2024 (doi: 10.1093/nar/gkae248)

# Rule 1: Run CMfinder to develop initial structure models from sequence alignments
rule cmfinder_initial:
    input:
        fasta="intergenic_regions/{gene_family}.fasta"
    output:
        sto="cmfinder_results/{gene_family}/initial/{gene_family}.sto",
        motif="cmfinder_results/{gene_family}/initial/{gene_family}.motif"
    params:
        outdir="cmfinder_results/{gene_family}/initial"
    threads: 4
    resources:
        mem_mb=8000,
        runtime=240  # 4 hours
    log:
        "logs/cmfinder/{gene_family}_initial.log"
    shell:
        """
        mkdir -p {params.outdir}
        
        docker run --rm \
            -v $(pwd)/{input.fasta}:/data/input.fasta:ro \
            -v $(pwd)/{params.outdir}:/data/output \
            cmfinder-container \
            cmfinder.pl \
            /data/input.fasta \
            > /data/output/{wildcards.gene_family}.sto \
            2> {log}
        
        # CMfinder also generates a motif file
        docker run --rm \
            -v $(pwd)/{params.outdir}:/data/output \
            cmfinder-container \
            sh -c "if [ -f /data/output/*.motif ]; then \
                mv /data/output/*.motif /data/output/{wildcards.gene_family}.motif; \
            fi"
        """

# Rule 2: Build covariance model using Infernal cmbuild
rule build_cm_model:
    input:
        sto="cmfinder_results/{gene_family}/initial/{gene_family}.sto"
    output:
        cm="cmfinder_results/{gene_family}/{gene_family}.cm"
    log:
        "logs/cmfinder/{gene_family}_cmbuild.log"
    shell:
        """
        docker run --rm \
            -v $(pwd)/{input.sto}:/data/input.sto:ro \
            -v $(pwd)/cmfinder_results/{wildcards.gene_family}:/data/output \
            cmfinder-container \
            cmbuild \
            /data/output/{wildcards.gene_family}.cm \
            /data/input.sto \
            2>&1 | tee {log}
        """

# Rule 3: Search for additional homologs using Infernal cmsearch
# This searches your database for more representatives
rule search_homologs:
    input:
        cm="cmfinder_results/{gene_family}/{gene_family}.cm",
        database="reference_genomes/refseq_bacteria.fasta"  # Your genomic database
    output:
        hits="cmfinder_results/{gene_family}/cmsearch/{gene_family}_hits.tblout",
        alignment="cmfinder_results/{gene_family}/cmsearch/{gene_family}_hits.sto"
    params:
        outdir="cmfinder_results/{gene_family}/cmsearch",
        evalue=0.01  # E-value threshold
    threads: 8
    resources:
        mem_mb=16000,
        runtime=480  # 8 hours
    log:
        "logs/cmfinder/{gene_family}_cmsearch.log"
    shell:
        """
        mkdir -p {params.outdir}
        
        docker run --rm \
            -v $(pwd)/{input.cm}:/data/model.cm:ro \
            -v $(pwd)/{input.database}:/data/database.fasta:ro \
            -v $(pwd)/{params.outdir}:/data/output \
            cmfinder-container \
            cmsearch \
            --cpu {threads} \
            -E {params.evalue} \
            --tblout /data/output/{wildcards.gene_family}_hits.tblout \
            -A /data/output/{wildcards.gene_family}_hits.sto \
            /data/model.cm \
            /data/database.fasta \
            2>&1 | tee {log}
        """

# Rule 4: Extract sequences from cmsearch hits
rule extract_hit_sequences:
    input:
        hits="cmfinder_results/{gene_family}/cmsearch/{gene_family}_hits.tblout",
        database="reference_genomes/refseq_bacteria.fasta"
    output:
        fasta="cmfinder_results/{gene_family}/expanded/{gene_family}_expanded.fasta"
    params:
        outdir="cmfinder_results/{gene_family}/expanded"
    log:
        "logs/cmfinder/{gene_family}_extract.log"
    shell:
        """
        mkdir -p {params.outdir}
        
        docker run --rm \
            -v $(pwd)/{input.hits}:/data/hits.tblout:ro \
            -v $(pwd)/{input.database}:/data/database.fasta:ro \
            -v $(pwd)/{params.outdir}:/data/output \
            cmfinder-container \
            esl-sfetch \
            --index /data/database.fasta \
            -f /data/hits.tblout \
            > /data/output/{wildcards.gene_family}_expanded.fasta \
            2> {log}
        """

# Rule 5: Refine alignment with CMfinder on expanded set
# This iteratively refines the structural model
rule cmfinder_refine:
    input:
        fasta="cmfinder_results/{gene_family}/expanded/{gene_family}_expanded.fasta"
    output:
        sto="cmfinder_results/{gene_family}/refined/{gene_family}_refined.sto",
        motif="cmfinder_results/{gene_family}/refined/{gene_family}_refined.motif"
    params:
        outdir="cmfinder_results/{gene_family}/refined"
    threads: 4
    resources:
        mem_mb=12000,
        runtime=360  # 6 hours
    log:
        "logs/cmfinder/{gene_family}_refine.log"
    shell:
        """
        mkdir -p {params.outdir}
        
        docker run --rm \
            -v $(pwd)/{input.fasta}:/data/input.fasta:ro \
            -v $(pwd)/{params.outdir}:/data/output \
            cmfinder-container \
            cmfinder.pl \
            /data/input.fasta \
            > /data/output/{wildcards.gene_family}_refined.sto \
            2> {log}
        
        docker run --rm \
            -v $(pwd)/{params.outdir}:/data/output \
            cmfinder-container \
            sh -c "if [ -f /data/output/*.motif ]; then \
                mv /data/output/*.motif /data/output/{wildcards.gene_family}_refined.motif; \
            fi"
        """

# Rule 6: Build final covariance model
rule build_final_cm:
    input:
        sto="cmfinder_results/{gene_family}/refined/{gene_family}_refined.sto"
    output:
        cm="cmfinder_results/{gene_family}/final/{gene_family}_final.cm"
    params:
        outdir="cmfinder_results/{gene_family}/final"
    log:
        "logs/cmfinder/{gene_family}_final_cm.log"
    shell:
        """
        mkdir -p {params.outdir}
        
        docker run --rm \
            -v $(pwd)/{input.sto}:/data/input.sto:ro \
            -v $(pwd)/{params.outdir}:/data/output \
            cmfinder-container \
            cmbuild \
            /data/output/{wildcards.gene_family}_final.cm \
            /data/input.sto \
            2>&1 | tee {log}
        """

# Rule 7: Calibrate the covariance model for accurate E-values
rule calibrate_cm:
    input:
        cm="cmfinder_results/{gene_family}/final/{gene_family}_final.cm"
    output:
        cm_calibrated="cmfinder_results/{gene_family}/final/{gene_family}_final_calibrated.cm"
    threads: 4
    resources:
        mem_mb=8000,
        runtime=120
    log:
        "logs/cmfinder/{gene_family}_calibrate.log"
    shell:
        """
        docker run --rm \
            -v $(pwd)/{input.cm}:/data/input.cm:ro \
            -v $(pwd)/cmfinder_results/{wildcards.gene_family}/final:/data/output \
            cmfinder-container \
            cmcalibrate \
            --cpu {threads} \
            /data/input.cm \
            2>&1 | tee {log}
        
        # Move calibrated model
        docker run --rm \
            -v $(pwd)/cmfinder_results/{wildcards.gene_family}/final:/data/output \
            cmfinder-container \
            cp /data/input.cm /data/output/{wildcards.gene_family}_final_calibrated.cm
        """

# Rule 8: Generate visualization with R2R (consensus structure)
rule visualize_structure:
    input:
        sto="cmfinder_results/{gene_family}/refined/{gene_family}_refined.sto"
    output:
        pdf="cmfinder_results/{gene_family}/visualization/{gene_family}_structure.pdf",
        svg="cmfinder_results/{gene_family}/visualization/{gene_family}_structure.svg"
    params:
        outdir="cmfinder_results/{gene_family}/visualization"
    log:
        "logs/cmfinder/{gene_family}_r2r.log"
    shell:
        """
        mkdir -p {params.outdir}
        
        docker run --rm \
            -v $(pwd)/{input.sto}:/data/input.sto:ro \
            -v $(pwd)/{params.outdir}:/data/output \
            cmfinder-container \
            r2r \
            --GSC-weighted-consensus \
            /data/input.sto \
            /data/output/{wildcards.gene_family}_structure.pdf \
            2>&1 | tee {log}
        
        # Also generate SVG
        docker run --rm \
            -v $(pwd)/{input.sto}:/data/input.sto:ro \
            -v $(pwd)/{params.outdir}:/data/output \
            cmfinder-container \
            r2r \
            --GSC-weighted-consensus \
            /data/input.sto \
            /data/output/{wildcards.gene_family}_structure.svg \
            2>&1 >> {log}
        """

# Master rule to run complete pipeline for all gene families
rule all:
    input:
        expand("cmfinder_results/{gene_family}/final/{gene_family}_final_calibrated.cm",
               gene_family=GENE_FAMILIES),
        expand("cmfinder_results/{gene_family}/visualization/{gene_family}_structure.pdf",
               gene_family=GENE_FAMILIES)

# Optional: Simple single-step CMfinder run (if you just want basic analysis)
rule cmfinder_simple:
    input:
        fasta="intergenic_regions/{gene_family}.fasta"
    output:
        sto="cmfinder_results_simple/{gene_family}.sto"
    threads: 4
    resources:
        mem_mb=8000
    log:
        "logs/cmfinder_simple/{gene_family}.log"
    shell:
        """
        mkdir -p cmfinder_results_simple
        
        docker run --rm \
            -v $(pwd)/{input.fasta}:/data/input.fasta:ro \
            -v $(pwd)/cmfinder_results_simple:/data/output \
            cmfinder-container \
            cmfinder.pl \
            /data/input.fasta \
            > /data/output/{wildcards.gene_family}.sto \
            2> {log}
        """