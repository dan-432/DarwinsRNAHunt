"""
author: Daniel Dachs
date: 11/02/2026
version: 1

Script with useful functions for querying interpro api.

Park, Y. M., Squizzato, S., Buso, N., Gur, T., & Lopez, R. (2017). The EBI search engine: 
EBI search as a service-making biological data accessible for all. Nucleic acids research, 
45(W1), W545-W549. https://doi.org/10.1093/nar/gkx359
"""

# standard library modules
import json, ssl
from urllib import request
from urllib.error import HTTPError
from time import sleep
from ete3 import NCBITaxa

#BASE_URL = "https://www.ebi.ac.uk:443/interpro/api/taxonomy/uniprot/protein/entry/InterPro/{prot_id}/?page_size=200"

def get_taxonomic_list(ip_protein_id: str, base_url: str):
    #disable SSL verification to avoid config issues
    context = ssl._create_unverified_context()

    # we need to iterate through sites of page size 200 (or whatever we've set)

    next = base_url.format(prot_id = ip_protein_id)

    results = []

    attempts = 0

    # until we get to the end of our results
    while next:
        try:
            req = request.Request(next, headers={"Accept": "application/json"})
            res = request.urlopen(req, context=context)

            # If the API times out due a long running query, this is best practice set out by interpro
            if res.status == 408:
                # wait just over a minute
                sleep(61)
                # then continue this loop with the same URL
                continue
            elif res.status == 204:
                #no data so leave loop
                break
            
            payload = json.loads(res.read().decode())

            if next != payload["next"]:
                next = payload["next"]

            # if we succeeded reset the attempts count
            attempts = 0

            results.extend(payload["results"])

        except HTTPError as e:
            # If there is a different HTTP error, it wil re-try 3 times before failing
            if attempts < 3:
                attempts += 1
                sleep(10)
                continue
            else:
                print("LAST URL: " + next)
                raise e

    return results


def download_and_save_taxanomic_info(ip_protein_id, file, base_url):
    tax_list = get_taxonomic_list(ip_protein_id, base_url)

    with open(file, "w") as f:
        json.dump(tax_list, f, indent=4)
    
    return tax_list

def load_taxonomic_info(file, children=False):
    """
    Load taxonomic info from a JSON file and return a list of taxa IDs.
    Args:
        file (str): Path to the JSON file containing taxonomic info.
        children (bool): Whether to include child taxa IDs.
    Returns:
        list: List of taxa IDs extracted from the JSON file."""
    info = []
    with open(file) as f:
        info = json.load(f)

    return get_taxa_ids(info, children=children)

def get_taxa_ids(taxa_info, children=False):
    taxa_ids = []

    for t in taxa_info:
        taxa_data = t["metadata"]
        taxa_ids.append(taxa_data["accession"])

        if children and taxa_data["children"]:
            # if we want to dig down to taxa children and these exist
            taxa_ids.extend(taxa_data["children"])


    return taxa_ids

def main():
    #download_and_save_taxanomic_info("IPR007024")
    taxa_info = load_taxonomic_info()

    taxa_ids = get_taxa_ids(taxa_info)

    print(len(taxa_ids))

    print("\nwith %d species" % len(get_taxa_ids(taxa_info, True)))

if __name__ == "__main__":
    main()