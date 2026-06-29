"""
author: Danil Dachs
date: 16/05/2026
version: 1

This module provides utilities for downloading and loading Rfam RNA family data, including covariance models (CMs) and Stockholm format alignments."""

from DarwinsRNAHunt.downloadallfastafiles import download_file
import gzip
import os

def download_gzipped_file(url, download_location):
    """
    Download gzip file to specified directory.
    
    Args:
        url (str): URL to the gzip file (e.g., Rfam.cm.gz)
        download_location (str): Local directory to save the downloaded file

    Returns:
        str: Local file path where the gzip file was saved
    """
    gzip_file = download_file(url, download_location+".temp.gz")
    with gzip.open(gzip_file, 'rb') as f_in:
        with open(download_location, 'wb') as f_out:
            f_out.write(f_in.read())

    os.remove(gzip_file)
    return download_location