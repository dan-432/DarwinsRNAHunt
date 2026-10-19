# DarwinsRNAHunt

A Snakemake-based honours project workflow for discovering RNA structures that are evolutionarily co-conserved in non-coding flanking regions of genes.

## Overview

This repository defines a pipeline that:

- downloads or consumes reference data from Rfam, GTDB, and InterPro;
- selects target genomes and taxonomic subsets;
- searches genome proteomes for protein domains of interest;
- extracts upstream and downstream flanking regions around domain hits;
- performs motif discovery and clustering for candidate RNA structures;
- includes dedicated control workflows for Rfam families and RMfam models.

## Pipeline flow

The main workflow entrypoint is `workflow/Snakefile`.
It reads parameters from `config/config.yaml` and includes rule definitions from `workflow/rules/`.

Typical pipeline stages:

1. `download_interpro_taxonomy` - fetch taxonomy mappings for the target InterPro IDs.
2. `download_gtdb_data` - download GTDB bacterial tree and metadata.
3. `download_rfam` - download Rfam covariance models and genome index data.
4. `select_genomes` - filter genomes by taxonomy, completeness, and contamination thresholds.
5. `download_genomes` - download selected genomes from NCBI into the local data store.
6. `locate_domains` - search genome proteomes for the target domain HMM profiles.
7. `extract_flanking_sequences` - extract 5' and 3' flanking sequences and CDS regions around hits.
8. `cmfinder` / `rfam_positive_controls` rules - run motif discovery and control analyses using additional rule modules.

The workflow is organized so that heavy data preparation and control workflows are split into `workflow/rules/*.smk`, while specialized Python and shell helper scripts live in `workflow/src/`.

## How to run

Snakemake is required, install following instructions from [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html). 
I recommend the Conda installation method.
envvionment.yaml in this directory defines a base conda environment for snakemake runs.

There is a wrapper script for Slurm-based execution:

```bash
./run_snakemake.sh
```

That script activates a Conda environment and runs Snakemake with:

- `--snakefile workflow/Snakefile`
- `--use-conda`
- `--conda-base-path ~/miniforge3`
- `--latency-wait 60`
- `--retries 2 --keep-going`

To run Snakemake directly from this repository:

```bash
cd /path/to/DarwinsRNAHunt
snakemake --snakefile workflow/Snakefile --cores 16 --use-conda --latency-wait 60 --retries 2 --keep-going
```

If Snakemake becomes locked, run `snakemake --unlock` before rerunning.

## Configuration

The main configuration files are `config/preferences.yaml`, `config/resources.yaml` and `config/technical.yaml`.

`config/preferences.yaml` is the start place for changing how the pipeline runs, it controls:

- target InterPro IDs, to set the targets for the pipeline, under `interpro.target_interpro_ids`
- HMMER search thresholds under `hmmer:`
- flanking region extraction settings under `flanking:`
- MMseqs2 clustering parameters under `clustering.mmseqs2:`
- control workflow settings under `controls:`

For resource control (changing the run time or memory alocation for particular steps) use `config/resources.yaml`.

Technical setting like web API URLs are set in `config/technical.yaml`.

## File structure

Top-level layout:

- `config/`
  - `config/preferences.yaml` - main Snakemake config file.
  - `config/resources.yaml`
  - `config/technical.yaml`
- `workflow/`
  - `Snakefile` - main workflow entrypoint.
  - `rules/` - Snakemake rule modules for pipeline stages.
  - `src/` - helper scripts used by rules.
  - `envs/` - Conda environment YAML files for rule execution.
- `results/` - generated workflow outputs.
- `logs/` - Snakemake and rule-level log files.
- `notebooks/` - analysis and debugging scripts and notebooks
- `data/` - reference data downloads and local asset files used by rules.
  - `data/genomes` - stores all genome files under accession folders, download checks for already existing genomes
- `run_snakemake.sh` - example wrapper script for cluster execution.
- `.tests` - unit tests

## Notes

- `results/` is the primary destination for downstream analysis outputs. Look at `results/03_manual_motif_discovery/**domain**_domain_flanking_**position**/final_report.html` for a final report of discovered candidate structured RNAs.
