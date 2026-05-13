"""
author: Daniel Dachs
date: 05/02/2026
version: 1

Script with useful functions for querying and retreiving genome sequences from the EBI API.

Park, Y. M., Squizzato, S., Buso, N., Gur, T., & Lopez, R. (2017). The EBI search engine: 
EBI search as a service-making biological data accessible for all. Nucleic acids research, 
45(W1), W545-W549. https://doi.org/10.1093/nar/gkx359
"""

import requests
import sys
import json

# Define the request URL
REQUEST_URL = "https://www.ebi.ac.uk/ena/portal/api/search"

# find taxa ID on ncbi
TAXAID = 1117 # cyanobacteria

# genome metadata fields of interest
FIELDS = "base_count,accession,genome_representation,tax_id,strain,scientific_name,last_updated,assembly_level,sample_accession,study_description"

def get_genome_len(genome):
    """Returns genome length from genome meta data.
    Parameters:
    genome (dict): genome meta data returned by EBI search of type json. feilds/keys expected is base_count
    
    Returns:
    int: genome base count"""
    
    return int(genome["base_count"])

def clean_genome_data(genome_data):
    """Filters genomes from EBI API for clean single strain assemblies. Note this is quite conservative. Bit of a yucky text based search but apparently not many other options:(
    Parameters:
    gemome (list[dict]) : list of genome meta data returned by EBI search of type json. feilds/keys expected are base_count,accession,genome_representation,tax_id,strain,scientific_name,last_updated,assembly_level,sample_accession,study_description
    
    Returns:
    list[dict]: cleaned list of genomes"""

    # study types we don't like (key words or partials)
    items_to_remove = ["complexes", "specimens", "strains", "collection","isolates","metagenom"]

    # keeps only genomes where the study description does not include any of our key words
    filtered_genomes =[g for g in genome_data if not any(ex in g["study_description"].lower() for ex in items_to_remove)]

    # sort by genome size, largest->smallest
    filtered_genomes.sort(key=get_genome_len, reverse=True)

    return filtered_genomes

def retreive_genomes(genome_ids, bacterial_genome_only):
    """TO BE CREATED"""

    return None

def retreive_genome_metadata(taxa_id, fields):

    retreive_metadata_payload = {
    "result": "assembly",
    "query": 'tax_tree(%d) AND genome_representation="full"' % (taxa_id),
    "limit": "0",
    "fields": fields,
    "format": "json"
    }

    headers = {
        "Content-Type": "application/x-www-form-urlencoded"
    }

    return requests.post(REQUEST_URL, data=retreive_metadata_payload, headers=headers).json()

def get_taxa_ids(genome_meta_data):

    taxa_ids = []

    for g in genome_meta_data:
        taxa_ids.append(g["tax_id"])

    return taxa_ids

def __main__():

    g = retreive_genome_metadata(TAXAID, FIELDS)

    genome_meta_data = clean_genome_data(g)

    print("%d genomes returned from taxa id %d" % (len(genome_meta_data), TAXAID))

    with open("bacteria-genomes.json", "w") as f:
        json.dump(genome_meta_data, f, indent=4)


if __name__ == "__main__":
    __main__()