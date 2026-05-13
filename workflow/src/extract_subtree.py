#!/usr/bin/env python3
"""
Extract subtree of domain containing species of specified taxa and save accessions
"""

import argparse
from pathlib import Path
from DarwinsRNAHunt.gtdb_access import load_tree, load_metadata, save_tree_image, get_ncbi_to_gtdb_dict, save_tree, save_metadata, save_genome_list
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
    
    args = parser.parse_args()
    
    # Load tree
    print(f"Loading tree from {args.input_tree}...")
    tree = load_tree(str(args.input_tree))

    print(f"Loading metadata from {args.metadata}...")
    metadata = load_metadata(str(args.metadata))

    if args.to_keep:
        # load taxanomic info from file and extract taxa ids
        ncbi_taxa = get_taxa_ids(load_taxonomic_info(str(args.to_keep)))

        print(f"Trimming tree to {len(ncbi_taxa)} NCBI taxa")

        # translate from ncbi to gtdb

        taxa_dict = get_ncbi_to_gtdb_dict(args.metadata)

        gtdb_taxaids = set()
        no_mapping_count = 0

        for ncbi_id in ncbi_taxa:
            if ncbi_id in taxa_dict.keys():
                gtdb_taxaids.update(taxa_dict[ncbi_id])
            else:
                no_mapping_count += 1

        print(f"NCBI taxonomy IDs which do not map to GTDB accessions: {no_mapping_count} of a total {len(ncbi_taxa)}")

        # trim tree
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
