#!/bin/bash
#SBATCH --partition=aoraki
#SBATCH --job-name=snakemake_dan
#SBATCH --output=snakemake_%j.log
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32          # Request more cores
#SBATCH --mem=128GB                  # Request more memory
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=dacda462@student.otago.ac.nz

source ~/miniforge3/etc/profile.d/conda.sh
export PYTHONNOUSERSITE=1
conda activate snakemake

# Run Snakemake locally on the compute node
snakemake \
  --snakefile workflow/Snakefile \
  --cores 32 \
  --use-conda \
  --conda-base-path ~/miniforge3 \
  --latency-wait 60 \
  --retries 2 --keep-going \
  --rerun-incomplete \
  results/00_controls/rfam/RF01315/07_cluster_analysis/phylo_report.tsv 

# --unlock \ add this line if you need to unlock the working directory