import argparse
from DarwinsRNAHunt.rfam_access import download_gzipped_file

def main():
    parser = argparse.ArgumentParser(
        description='Download Rfam file'
    )
    parser.add_argument('--rfam-cm-url', required=True,
                       help='URL to Rfam file (e.g., Rfam.cm.gz)')
    parser.add_argument('--output-cm', required=True,
                       help='Output file path for downloaded file')
    parser.add_argument('--rfam-genome-index-url', required=True,
                       help='URL to Rfam genome index file (e.g., Rfam.full_region.gz)')
    parser.add_argument('--output-genome-index', required=True,
                       help='Output file path for downloaded genome index file')

    args = parser.parse_args()

    print(f"Downloading Rfam file from {args.rfam_cm_url} to {args.output_cm}")
    print(f"Downloading Rfam genome index file from {args.rfam_genome_index_url} to {args.output_genome_index}")

    download_gzipped_file(args.rfam_cm_url, args.output_cm)
    download_gzipped_file(args.rfam_genome_index_url, args.output_genome_index)


if __name__ == "__main__":
    main()