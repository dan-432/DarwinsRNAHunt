# Configuration parameters
MULTI_CM_FILE = "resources/rmfam/RMfam.cm"  # Path to the multi-CM file from Rfam
NUM_SEQS_TO_KEEP = 50          
INITIAL_POOL_SIZE = 500        
BIT_SCORE_THRESHOLD = 20.0     


# 1. This rule runs inside Conda to safely extract model names
checkpoint list_models:
    input:
        multi_cm = MULTI_CM_FILE
    output:
        list_txt = "results/00_controls/rmfam/model_list.txt"
    conda:
        "../envs/rna_motif_env.yaml"  # Path to your Conda env file containing infernal
    shell:
        """
        cmstat {input.multi_cm} | grep -v '^#' | awk 'NF>1 {{print $2}}' > {output.list_txt}
        """

rule extract_single_cm:
    input:
        multi_cm = MULTI_CM_FILE
    output:
        single_cm = temp("results/00_controls/rmfam/{model}.cm")
    conda:
        "../envs/rna_motif_env.yaml"
    shell:
        "cmfetch {input.multi_cm} {wildcards.model} > {output.single_cm}"

rule generate_and_filter_sequences:
    input:
        cm = "results/00_controls/rmfam/{model}.cm"
    output:
        pool = temp("results/00_controls/rmfam/{model}_raw_pool.fa"),
        tbl = temp("results/00_controls/rmfam/{model}_search_results.tbl"),
        filtered = "results/00_controls/rmfam/{model}_filtered.fa"
    params:
        n_keep = NUM_SEQS_TO_KEEP,
        pool_size = INITIAL_POOL_SIZE,
        threshold = BIT_SCORE_THRESHOLD
    conda:
        "../envs/rna_motif_env.yaml"
    shell:
        """
        cmemit -N {params.pool_size} {input.cm} > {output.pool}
        cmsearch --tblout {output.tbl} {input.cm} {output.pool} > /dev/null

        python3 -c "
import sys
passed_seqs = set()
with open('{output.tbl}', 'r') as f:
    for line in f:
        if line.startswith('#'): continue
        parts = line.split()
        if len(parts) > 14 and float(parts[14]) >= {params.threshold}:
            passed_seqs.add(parts[0])

kept_count = 0
with open('{output.pool}', 'r') as infile, open('{output.filtered}', 'w') as outfile:
    current_seq_name = ''
    current_seq_lines = []
    for line in infile:
        if line.startswith('>'):
            if current_seq_name in passed_seqs and kept_count < {params.n_keep}:
                outfile.write(''.join(current_seq_lines))
                kept_count += 1
            current_seq_name = line.strip().split()[0][1:]
            current_seq_lines = [line]
        else:
            current_seq_lines.append(line)
    if current_seq_name in passed_seqs and kept_count < {params.n_keep}:
        outfile.write(''.join(current_seq_lines))
        kept_count += 1
"
        """

# 2. This function triggers after the checkpoint finishes, reading the models safely
def aggregate_models(wildcards):
    checkpoint_output = checkpoints.list_models.get(**wildcards).output.list_txt
    with open(checkpoint_output) as f:
        models = [line.strip() for line in f if line.strip()]
    return expand("results/00_controls/rmfam/{model}_filtered.fa", model=models)

# 3. Final consumer rule that relies on the dynamic list of outputs
rule combine_all:
    input:
        aggregate_models
    output:
        "results/00_controls/rmfam_control_sequences.fa"
    shell:
        "cat {input} > {output}"
