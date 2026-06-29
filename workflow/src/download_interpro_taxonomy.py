#!/usr/bin/env python3
"""
Download taxonomic information from InterPro for a given domain
"""

import argparse
from DarwinsRNAHunt.interpro_access import download_and_save_taxanomic_info


def main():
    print("Starting up download_interpro_taxonomy.py...")
    parser = argparse.ArgumentParser(
        description='Download InterPro taxonomy data'
    )
    parser.add_argument('--interpro-id', required=True,
                       help='InterPro ID (e.g., IPR007024)')
    parser.add_argument('--interpro-tax-base-url', required=True,
                       help='Base url for interpro, e.g. https://www.ebi.ac.uk:443/interpro/api/taxonomy/uniprot/protein/entry/InterPro/{prot_id}/?page_size=200')
    parser.add_argument('--output', required=True,
                       help='Output JSON file')
    
    args = parser.parse_args()
    
    # Download taxonomy data
    print(f"Downloading taxonomy for {args.interpro_id}...")
    download_and_save_taxanomic_info(args.interpro_id, args.output, args.interpro_tax_base_url)
    
    print(f"Saved taxonomy data to {args.output}")


if __name__ == '__main__':
    main()
