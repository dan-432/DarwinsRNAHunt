#!/usr/bin/env python3
"""
Download GTDB tree and metadata - doing this in one script as they must be complementary
"""

import argparse
from DarwinsRNAHunt.gtdb_access import download_gtdb_bacteria_tree, download_metadata


def main():
    parser = argparse.ArgumentParser(
        description='Download tree and associated metadata from GTDB'
    )
    parser.add_argument('--tree-url', required=True,
                       help='URL to GTDB tree')
    parser.add_argument('--metadata-url', required=True,
                       help='URL to GTDB metadata - MUST MATCH TREE')
    parser.add_argument('--tree-output', required=True,
                       help='Output taxanomic tree file')
    parser.add_argument('--metadata-output', required=True,
                       help='Output taxanomic metadata file')
    
    args = parser.parse_args()

    print(f"Downloading taxanomic tree from {args.tree_url} to {args.tree_output}")
    
    download_gtdb_bacteria_tree(args.tree_url, args.tree_output)

    print(f"Downloading taxanomic tree metadata from {args.metadata_url} to {args.metadata_output}")
    
    download_metadata(args.metadata_url, args.metadata_output)


if __name__ == '__main__':
    main()
