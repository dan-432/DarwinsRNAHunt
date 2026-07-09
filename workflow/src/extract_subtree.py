#!/usr/bin/env python3
"""
Extract subtree of domain containing species of specified taxa and save accessions
"""

import argparse
import json
from pathlib import Path
from DarwinsRNAHunt.gtdb_access import load_tree, load_metadata, save_tree_image, get_ncbi_to_gtdb_dict, save_tree, save_metadata, save_genome_list, get_assembly_to_gtdb_dict, normalize_accession
from DarwinsRNAHunt.gtdb_taxonomy_analysis import trim_tree_to_taxa, trim_tree_to_order, get_tree_accession_ids
from DarwinsRNAHunt.interpro_access import load_taxonomic_info, get_taxa_ids 

def main():
    parser = argparse.ArgumentParser(
        description='Extract taxonomic subtree, note no filtering will occur unless optional arguments are defined.'
    )
    parser.add_argument('--input-tree', required=True, type=Path,
                       help='Input taxanomic tree file, Newick format')
    parser.add_argument('--metadata', required=True, type=Path,
                       help='GTDB metadata file, gzip csv')
    parser.add_argument('--to-keep', required=False, type=Path,
                        help='List of NCBI genome assemblies to keep, e.g. all of those containing domain of interest from Interpro')
    parser.add_argument('--target-taxanomic-order', required=False,
                       help='Target taxanomic order name GTDB nomenclature (e.g., 36.0:f__Burkholderiaceae)')
    parser.add_argument('--completeness', type=float, required=False,
                       help='Minimum completeness threshold')
    parser.add_argument('--contamination', type=float, required=False,
                       help='Maximum contamination threshold')
    parser.add_argument('--output-tree', required=True,
                        help='Output tree file')
    parser.add_argument('--output-metadata', required=True,
                       help='Output gzip CSV with filtered accession metadata')
    parser.add_argument('--output-genomes', required=True,
                       help='Output text file with list of filtered genome accession IDs', )
    parser.add_argument('--output-image', required=False,
                       help='Output tree image')
    parser.add_argument('--from-interpro', action='store_true',
                       help='If set, will load taxanomic info from Interpro JSON file')
    
    args = parser.parse_args()
    
    # Load tree
    print(f"Loading tree from {args.input_tree}...")
    tree = load_tree(str(args.input_tree))

    print(f"Loading metadata from {args.metadata}...")
    metadata = load_metadata(str(args.metadata))

    if args.to_keep:
        print(f"Loading taxa to keep from {args.to_keep}...")

        if args.from_interpro:
            print("InterPro input detected; extracting NCBI taxa IDs...")
            taxa_to_keep = load_taxonomic_info(str(args.to_keep))
            accession_mapping = get_ncbi_to_gtdb_dict(args.metadata)

            print(f"Trimming tree to {len(taxa_to_keep)} target taxa")

            gtdb_taxaids = set()
            unmapped_count = 0

            for taxon_id in taxa_to_keep:
                if taxon_id in accession_mapping:
                    gtdb_taxaids.update(accession_mapping[taxon_id])
                else:
                    unmapped_count += 1

        if not args.from_interpro:
            taxa_to_keep = []
            # if not from interpro, then we assume the user has provided a list of assembly accession IDs to keep
            with open(args.to_keep) as f:
                taxa_to_keep = json.load(f)

            print("Using assembly accession IDs...")
            accession_mapping = get_assembly_to_gtdb_dict(args.metadata)

            print(f"Trimming tree to {len(taxa_to_keep)} target taxa")

            gtdb_taxaids = set()
            unmapped_count = 0

            for taxon_id in taxa_to_keep:
                normalized_id = normalize_accession(taxon_id)
                if normalized_id in accession_mapping:
                    gtdb_taxaids.update(accession_mapping[normalized_id])
                else:
                    unmapped_count += 1

        print(
            f"{unmapped_count} of {len(taxa_to_keep)} IDs could not be mapped to GTDB accessions"
        )
        tree = trim_tree_to_taxa(tree, gtdb_taxaids)

    if args.target_taxanomic_order:
        # Find target node, i.e. narrow dowm to 
        print(f"Searching for node: {args.target_taxanomic_order}")
        tree = trim_tree_to_order(tree, args.target_taxanomic_order)
        
        if tree is None:
            raise ValueError(f"Node {args.target_taxanomic_order} not found in tree")
        
    # finished now with tree manipulations

    # save tree
    print(f"Saving tree to {args.output_tree}")

    save_tree(tree, args.output_tree)

    # Render tree
    # if args.output_image:
    #     print(f"Rendering tree to {args.output_image}...")
    #     save_tree_image(tree, args.output_image)
    #     print("Done")

        
    accession_names = get_tree_accession_ids(tree)
    
    # filter genomes to specified contamination and completeness

    select_accessions_metadata = metadata[metadata["accession"].isin(accession_names)]

    if args.completeness and args.contamination:
        select_accessions_metadata = select_accessions_metadata[(select_accessions_metadata["checkm2_completeness"]>args.completeness) & (select_accessions_metadata["checkm2_contamination"]<args.contamination)]
    elif args.completeness:
        select_accessions_metadata = select_accessions_metadata[(select_accessions_metadata["checkm2_completeness"]>args.completeness)]
    elif args.contamination:
        select_accessions_metadata = select_accessions_metadata[(select_accessions_metadata["checkm2_contamination"]<args.contamination)]


    print(f"Saving filtered genome metadata to {args.output_metadata}...")

    save_metadata(select_accessions_metadata, args.output_metadata)

    print(f"Saving list of {len(select_accessions_metadata)}genome acccession IDs to {args.output_genomes}...")

    save_genome_list(select_accessions_metadata, args.output_genomes)

    print("DONE")


if __name__ == '__main__':
    main()
