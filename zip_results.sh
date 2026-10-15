#!/bin/bash
#SBATCH --partition=aoraki
#SBATCH --job-name=snakemake_dan
#SBATCH --output=snakemake_%j.log
#SBATCH --error=snakemake_%j.err
#SBATCH --time=12:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8GB                  # Request more memory
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=dacda462@student.otago.ac.nz

tar -czvf discovery_results.tar.gz results/ ; tar -czvf discovery_data.tar.gz data/