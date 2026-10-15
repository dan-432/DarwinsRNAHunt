#!/usr/bin/env python3
"""Search genome FASTA files with an Infernal covariance model."""

import argparse
import json
import shutil
import subprocess
import tempfile
from pathlib import Path


def parse_args():
	parser = argparse.ArgumentParser(
		description="Run cmsearch across genomes listed in a metadata JSON file."
	)
	parser.add_argument(
		"metadata",
		type=Path,
		help="JSON file containing an 'accessions' list",
	)
	parser.add_argument("cm", type=Path, help="Covariance model to search with")
	parser.add_argument(
		"--genome-root",
		type=Path,
		default=Path("data/genomes"),
		help="Directory containing <accession>/genomic.fna (default: data/genomes)",
	)
	parser.add_argument(
		"--threads",
		type=int,
		default=1,
		help="Number of CPU threads passed to cmsearch (default: 1)",
	)
	parser.add_argument(
		"--evalue",
		type=float,
		help="Optional reporting E-value threshold passed to cmsearch with -E",
	)
	return parser.parse_args()


def load_accessions(metadata_path):
	with metadata_path.open() as metadata_file:
		metadata = json.load(metadata_file)

	accessions = metadata.get("accessions")
	if not isinstance(accessions, list) or not all(
		isinstance(accession, str) and accession for accession in accessions
	):
		raise ValueError("metadata must contain a non-empty-string 'accessions' list")
	return accessions


def count_tblout_hits(tblout_path):
	with tblout_path.open() as tblout_file:
		return sum(
			1
			for line in tblout_file
			if line.strip() and not line.lstrip().startswith("#")
		)


def search_genome(cm_path, fasta_path, threads, evalue, tblout_path):
	command = [
		"cmsearch",
		"--cpu",
		str(threads),
		"--noali",
		"--nohmmonly",
		"--tblout",
		str(tblout_path),
	]
	if evalue is not None:
		command.extend(["-E", str(evalue)])
	command.extend([str(cm_path), str(fasta_path)])
	subprocess.run(command, check=True)
	return count_tblout_hits(tblout_path)


def main():
	args = parse_args()
	if args.threads < 1:
		raise ValueError("--threads must be at least 1")
	if shutil.which("cmsearch") is None:
		raise FileNotFoundError("cmsearch was not found in PATH")
	if not args.metadata.is_file():
		raise FileNotFoundError(f"metadata file not found: {args.metadata}")
	if not args.cm.is_file():
		raise FileNotFoundError(f"covariance model not found: {args.cm}")

	accessions = load_accessions(args.metadata)
	total_hits = 0
	with tempfile.TemporaryDirectory(prefix="cmsearch_") as temp_dir:
		for accession in accessions:
			fasta_path = args.genome_root / accession / "genomic.fna"
			if not fasta_path.is_file():
				raise FileNotFoundError(f"genome FASTA not found: {fasta_path}")
			tblout_path = Path(temp_dir) / f"{accession}.tblout"
			hit_count = search_genome(
				args.cm, fasta_path, args.threads, args.evalue, tblout_path
			)
			total_hits += hit_count
			print(f"{accession}\t{hit_count}")

	print(f"Total hits\t{total_hits}")


if __name__ == "__main__":
	try:
		main()
	except (FileNotFoundError, ValueError, json.JSONDecodeError, subprocess.CalledProcessError) as error:
		raise SystemExit(f"error: {error}") from error
