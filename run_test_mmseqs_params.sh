#!/bin/bash
#SBATCH --partition=aoraki
#SBATCH --job-name=test_mmseqs_params
#SBATCH --output=test_mmseqs_params_%j.log
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32GB

set -euo pipefail

source ~/miniforge3/etc/profile.d/conda.sh
export PYTHONNOUSERSITE=1

conda run -n domain-taxonomy python notebooks/test_mmseqs_params.py \
    --input results/02_sequence_selection/cluster_1341/03_combined/cluster_1341_domain_flanking_upstream_filtered.fasta \
    --cm data/rfam/covariance_models/RF00050.cm \
    --outdir notebooks/mmseqstestresults/trimmed/ \
    --min-cluster-size 5
