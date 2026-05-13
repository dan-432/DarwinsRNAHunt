""" 
author: Daniel Dachs
date: 10/03/2026

Program for scraping fasta sequence files froma. website 
"""


import errno
import gzip
import json
import os
import time
import requests
from Bio import SeqIO

# Target web page URL (modify as needed)
#BASE_URL = "https://riboswitch.ribocentre.org/downloads/sequences"
#https://riboswitch.ribocentre.org/downloads/sequences/RF00050.fa.gz


def download_file(url, output):
    """Save the file from the link.
    
    Args:
        url: String link to file
        output: String path/to/file"""

    file_name = url.split("/")[-1]

    print(f"Downloading: {url}")
    try:
        response = requests.get(url, stream=True, timeout=10)
        with open(output, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        print(f"Saved: {output}")
    except requests.RequestException as e:
        print(f"Error downloading {url}: {e}")

    return output


def read_fasta_gz_sequences(filename):
    """
    Reads sequences from a gzipped FASTA file using Biopython.

    Args:
        filename (str): The path to the .fa.gz file.

    Returns:
        A list of sequence records (see SeqIO)
    """
    sequences = []
    # Open the gzipped file in text mode ('rt')
    with gzip.open(filename, "rt") as handle:
        # Use SeqIO.parse to iterate over the records in the FASTA file
        for record in SeqIO.parse(handle, "fasta"):
            
            sequences.append(record)
            # You can also access other attributes like:
            #print(f"ID: {record.id}")
            #print(f"Description: {record.description}")
            #print(f"Sequence Length: {len(record.seq)}")
            
    return sequences

def save_fasta_to_text(filename, sequences=None):
    """ Generates readable txt file from fasta.gz
    
    Args:
        filename (str): path to .fa.gz file.
        sequences(list[SeqIO record]): list of sequences with info, if null will be generated from fasta
    
    Returns:
        path to newly generated txt file
    """

    txt_path = filename + ".txt"

    if not sequences:
        sequences = read_fasta_gz_sequences(filename)

    with open(txt_path, "w", encoding="utf-8") as f:
            for s in sequences:
                # 'record.seq' contains the sequence object
                seq_info = str(s.id) + ": " + s.description + "\n" + str(s.seq) + "\n"
                f.write(seq_info)

    return txt_path





def test():
    # was too much work to automate, check https://riboswitch.ribocentre.org/sequences/ if these are up to date
    download_link_ids = [
                            "RF00050",
                            "RF00059",
                            "RF00080",
                            "RF00162",
                            "RF00167",
                            "RF00168",
                            "RF00174",
                            "RF00230",
                            "RF00234",
                            "RF00379",
                            "RF00380",
                            "RF00442",
                            "RF00504",
                            "RF00521",
                            "RF00522",
                            "RF00634",
                            "RF01051",
                            "RF01054",
                            "RF01055",
                            "RF01056",
                            "RF01057",
                            "RF01068",
                            "RF01482",
                            "RF01689",
                            "RF01704",
                            "RF01725",
                            "RF01727",
                            "RF01734",
                            "RF01739",
                            "RF01750",
                            "RF01763",
                            "RF01764",
                            "RF01767",
                            "RF01786",
                            "RF01826",
                            "RF01831",
                            "RF02680",
                            "RF02683",
                            "RF02885",
                            "RF02974",
                            "RF02977",
                            "RF03013",
                            "RF03038",
                            "RF03054",
                            "RF03057",
                            "RF03071",
                            "RF03165",
                            "RF03167",
                            "RF03168",
                            "RF03169",
                            "RF03170"
                            ]#get_download_links(BASE_URL)#

    download_dir = "data/riboswitch_sequences/riboswitches."

    try:
        os.makedirs(download_dir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    riboswitch_meta_data = {}

    for id in download_link_ids:

        link = "https://riboswitch.ribocentre.org/downloads/sequences/%s.fa.gz" % id
        seq_data = []

        p = download_file(link, download_dir)

        seqs = read_fasta_gz_sequences(p)

        for s in seqs:
            d = {}
            d["id"] = s.id
            d["accession"] = s.id.split("/")[0]
            d["file"] = p
            d["annotations"] = s.annotations

            seq_data.append(d)

        print(save_fasta_to_text(p, seqs))

        riboswitch_meta_data[id] = seq_data

    with open(os.path.join(download_dir, "meta_data.json"), "w") as f:
        json.dump(riboswitch_meta_data, f, indent=4)

if __name__ == "__main__":
    test()