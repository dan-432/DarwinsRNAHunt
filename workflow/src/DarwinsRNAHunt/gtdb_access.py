"""
author: Daniel Dachs
date: 16/03/2026
version: 2

This module provides utilities for downloading and loading the GTDB (Genome Taxonomy Database)
bacterial phylogenetic tree and associated metadata. It includes helper functions for mapping
between NCBI taxonomy IDs and GTDB accessions.

Data sources:
- GTDB data, this uses latest: https://gtdb.ecogenomic.org/
"""

from DarwinsRNAHunt.downloadallfastafiles import download_file
from ete3 import Tree
import pandas
import json

# GTDB data URLs - latest release
#gtdb_bacteria_tree_url = "https://data.gtdb.ecogenomic.org/releases/latest/bac120.tree"
#gtdb_bacteria_metadata_url = "https://data.gtdb.ecogenomic.org/releases/latest/bac120_metadata.tsv.gz"

# Local storage configuration
#download_location = "data/taxa_trees"
#tree_file_name = "bac120.tree"
#metadata_file_name = "bac120_metadata.tsv.gz"


def download_gtdb_bacteria_tree(gtdb_bacteria_tree_url, download_location):
    """
    Download the GTDB bacterial phylogenetic tree (bac120.tree) to specified directory.
    
    The tree contains ~85,000+ bacterial genomes organized by phylogenetic relationships.
    
    Returns:
        str: Local file path where the tree was saved
    """
    return download_file(gtdb_bacteria_tree_url, download_location)

def save_tree(tree, file_name):
    """
    Save phylogenetic tree to data folder using ete3.
    
    Args:
        tree: ete3.Tree object
        file_name: Name of file to save

    Returns:
        str: Local file path where the tree was saved
    """
    # format=1 preserves branch lengths
    tree.write(outfile=file_name, format=1, format_root_node=True)
    return file_name

def load_tree(file_name):
    """
    Load a GTDB bacterial tree from a Newick format file using ete3.
    
    Args:
        file_name (str): Path to the .tree file (Newick format)
        
    Returns:
        ete3.Tree: Phylogenetic tree object
        
    Note:
        GTDB trees have quoted node names (e.g., 'RS_GCF_000006965.1').
        format=1 allows Newick with branch lengths.
        quoted_node_names=True handles GTDB's naming convention.
    """
    return Tree(file_name, format=1, quoted_node_names=True)


def download_and_load_bacteria_tree(gtdb_bacteria_tree_url, download_location):
    """
    Download and load the GTDB bacterial tree in one step.
    
    Convenience function that combines download and loading operations.
    
    Args:
        use_ete3 (bool): If True, use ete3 (recommended). If False, use Bio.Phylo.
    
    Returns:
        ete3.Tree or Bio.Phylo.Newick.Tree: Loaded phylogenetic tree
    """
    location = download_gtdb_bacteria_tree(gtdb_bacteria_tree_url, download_location)
    
    return load_tree(location)

def save_tree_image(tree, output_image):
    """
     Save tree as schematic image.
     
     Args:
        tree: ete3.Tree object
        output_image: string output image file, e.g. "results/trees/taxa.png"
    """
     
    # maybe add fomatting and change colours etc. not urgent
    tree.render(output_image)

def load_metadata(metadata_file):
    """
    Load GTDB bacterial metadata from a TSV file.
    
    Args:
        file_name (str): Path to the metadata file (.tsv.gz or .tsv)
        
    Returns:
        pandas.DataFrame: Metadata table with columns including:
            - accession: GTDB genome accession (e.g., GCA_948625425.1)
            - ncbi_taxid: NCBI taxonomy ID
            - gtdb_taxonomy: GTDB taxonomic classification string
            - checkm_completeness: Genome completeness percentage
            - checkm_contamination: Estimated contamination
            - and many more...
            
    Note:
        Pandas automatically handles .gz decompression when reading
    """
    dtype_mapping = {'accession': 'string', 'ncbi_taxid': 'string'}
    # Pandas read_csv automatically decompresses .gz files
    return pandas.read_csv(metadata_file, sep="\t", dtype=dtype_mapping)

def download_metadata(metadata_url, download_dir):
    """
    Download the GTDB bacterial metadata in one step.
    
    Convenience function that combines download and loading operations.
    
    Returns:
        String: path to file
    """
    return download_file(metadata_url, download_dir)


def download_and_load_metadata(metadata_url, download_dir):
    """
    Download and load the GTDB bacterial metadata in one step.
    
    Convenience function that combines download and loading operations.
    
    Returns:
        pandas.DataFrame: Loaded metadata table
        
    Example:
        >>> metadata = download_and_load_bacteria_metadata()
        >>> print(metadata.columns.tolist())
    """
    location = download_metadata(metadata_url, download_dir)
    return load_metadata(location)


def save_metadata(metadata_frame, metadata_file):
    """
    Save metadata to compressed tsv.gz file.

    Args: 
        metadata_frame (pandas.DataFrame): Metadata table with columns including:
            - accession: GTDB genome accession (e.g., GCA_948625425.1)
            - ncbi_taxid: NCBI taxonomy ID
            - gtdb_taxonomy: GTDB taxonomic classification string
            - checkm_completeness: Genome completeness percentage
            - checkm_contamination: Estimated contamination
            - and many more...
        metadata_file (String): path to save location, e.g. "data/taxa_trees/" or absolute
    """

    metadata_frame.to_csv(metadata_file, index=False, sep="\t", compression = 'gzip')

def save_genome_list(metadata_frame: pandas.DataFrame, output):
    """
    Save list NCBI genome accession IDs to text file from given metadata.

    Args: 
        metadata_frame (pandas.DataFrame): Metadata table with columns including:
            - accession: GTDB genome accession (e.g., GCA_948625425.1)
            - ncbi_taxid: NCBI taxonomy ID
            - gtdb_taxonomy: GTDB taxonomic classification string
            - checkm_completeness: Genome completeness percentage
            - checkm_contamination: Estimated contamination
            - and many more...
        output (String): path to save location, e.g. "data/taxa_trees/file.json" or absolute
    """

    try:

        acc_ids = metadata_frame["ncbi_genbank_assembly_accession"]

        acc_ids.to_json(output, orient='records', indent=4)

    except Exception as e:
        print("An unexpeccted error ocured while saving taxanomic meta data")
        raise e
    finally:
        print(f"Successfully saved to {output}")

def load_genome_list(file):
    """
    Save list NCBI genome accession IDs to text file from given metadata.

    Args: 
        file (String): path to json file, e.g. "data/taxa_trees/file.json" or absolute
    """
    try:
        with open(file, 'r') as json_file:
            data_list = json.load(json_file)
    except FileNotFoundError as e:
        # This block executes if the file is not found
        print(f"Error: The file '{file}' was not found.")
        raise e
    except PermissionError as e:
        # This block handles cases where you don't have permission to access the file
        print(f"Error: You do not have permission to open '{file}'.")
        raise e
    except Exception as e:
        # This general exception block catches any other potential errors
        print(f"An unexpected error occurred: {e}")
        raise e

    return data_list


# ===== Helper Functions =====

def get_ncbi_to_gtdb_dict(metadata_file, loadfile=True):
    """
    Create a mapping from NCBI taxonomy IDs to GTDB accession numbers.
    
    This is useful for converting between NCBI taxonomic identifiers (used by
    many databases like InterPro, UniProt) and GTDB accession numbers (used in
    the phylogenetic tree).
    
    Args:
        metadata_file (String): 
        loadfile (bool, optional): If True, load from existing local file. 
                              If False, download fresh metadata. Default: True
        
    Returns:
        dict: Dictionary mapping NCBI taxonomy IDs to lists of GTDB accessions
            Format: {ncbi_taxid: [accession1, accession2, ...]}
            
    Example:
        >>> ncbi_to_gtdb = get_ncbi_to_gtdb_dict()
        >>> ncbi_to_gtdb[562]  # E. coli tax ID
        ['GCA_000005845.2', 'GCA_000008865.2', ...]
        
    Note:
        One NCBI taxonomy ID can map to multiple GTDB accessions because GTDB reclassifies
        
    Raises:
        Exception: If loading fails, suggests trying with load=False to re-download
    """
    if loadfile:
        # Try to load from existing local file
        try:
            md = load_metadata(metadata_file)
        except Exception as e:
            print("TRY DOWNLOADING METADATA AGAIN WITH load=False PARAMETER, "
                  "ALSO DOUBLE CHECK DOWNLOAD LOCATION")
            raise e
    else:
        # Download fresh copy of metadata
        md = download_and_load_metadata(metadata_file)
    
    # Group accessions by NCBI taxonomy ID
    # Result: {taxid: [list of accessions]}
    taxa_dict = md.groupby('ncbi_taxid')['accession'].apply(list).to_dict()
    
    return taxa_dict

def normalize_accession(accession):
    """
    Normalize genome accession by stripping prefixes, GCF_, GCA_, RS_, GB_
    Args:
        accession (str): Genome accession string
    Returns:
        str: Normalized accession string
    """
    accession = str(accession).strip()
    prefixes = ("RS_", "GB_", "GCF_", "GCA_")
    changed = True
    while changed:
        changed = False
        for prefix in prefixes:
            if accession.startswith(prefix):
                accession = accession[len(prefix):]
                changed = True
    return accession

def get_assembly_to_gtdb_dict(metadata_file, loadfile=True):
    """
    Create a mapping from NCBI assembly accessions to representative GTDB accession numbers.
    
    This is useful for converting between NCBI assembly identifiers (used by
    many databases like InterPro, UniProt) and GTDB accession numbers (used in
    the phylogenetic tree).
    
    Args:
        metadata_file (String): 
        loadfile (bool, optional): If True, load from existing local file. 
                              If False, download fresh metadata. Default: True
        
    Returns:
        dict: Dictionary mapping NCBI assembly accessions to lists of GTDB accessions
            Format: {assembly_accession: [accession1, accession2, ...]}
            
    Example:
        >>> assembly_to_gtdb = get_assembly_to_gtdb_dict()
        >>> assembly_to_gtdb['GCA_000005845.2']  # E. coli assembly accession
        ['GCA_000005845.2', 'GCA_000008865.2', ...]
        
    Note:
        One NCBI assembly accession can map to multiple GTDB accessions because GTDB reclassifies
        
    Raises:
        Exception: If loading fails, suggests trying with load=False to re-download
    """
    if loadfile:
        try:
            md = load_metadata(metadata_file)
        except Exception as e:
            print("TRY DOWNLOADING METADATA AGAIN WITH load=False PARAMETER, "
                  "ALSO DOUBLE CHECK DOWNLOAD LOCATION")
            raise e
    else:
        md = download_and_load_metadata(metadata_file)

    md = md.dropna(subset=["ncbi_genbank_assembly_accession"]).copy()
    md["_normalized_key"] = md["ncbi_genbank_assembly_accession"].apply(normalize_accession)

    # keys: normalized assembly accession (matches GCA_ or GCF_ query IDs)
    # values: original GTDB 'accession' column, RS_/GB_ prefix intact —
    # needed as-is later since that's the format tree leaf names use
    assembly_dict = md.groupby("_normalized_key")["gtdb_genome_representative"].apply(list).to_dict()

    return assembly_dict

def ncbi_accessions_to_gtdb(accessions, assembly_to_gtdb_dict):
    """
    Map a list of NCBI assembly accessions to a flat, deduplicated list of
    GTDB representative accessions (RS_/GB_ prefix intact, tree-leaf format).

    Args:
        accessions: list of NCBI assembly accessions, e.g. ['GCA_000005845.2', ...]
        assembly_to_gtdb_dict: dict from your get_assembly_to_gtdb_dict()
            {normalized_ncbi_accession: [gtdb_genome_representative, ...]}

    Returns:
        (gtdb_accessions, unmapped): flat sorted list of unique GTDB accessions,
        and the list of input accessions that had no entry in the dict
    """
    gtdb_accessions = set()
    unmapped = []

    for accession in accessions:
        normalized = normalize_accession(accession)
        matches = assembly_to_gtdb_dict.get(normalized)
        if matches:
            gtdb_accessions.update(matches)
        else:
            unmapped.append(accession)

    if unmapped:
        print(f"{len(unmapped)}/{len(accessions)} accessions had no GTDB mapping "
              f"(showing up to 10): {unmapped[:10]}")

    return sorted(gtdb_accessions), unmapped
