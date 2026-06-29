#!/usr/bin/env python3
"""
Generate ESM2 embeddings for protein sequences (per-family).

Header format: >GCF_002897375.1|WP_103129364.1|CDS|NZ_BFAG01000006.1:207150-207878(+)

Usage:
    python embed_proteins_esm2.py \
        --fasta combined_flanking_cds.faa \
        --output embeddings.npz \
        --device cuda  # or cpu
"""

import argparse
import torch
import numpy as np
import pandas as pd
from pathlib import Path
from tqdm import tqdm
import re

try:
    import esm
except ImportError:
    print("Error: ESM not installed. Install with:")
    print("  pip install fair-esm2")
    exit(1)

def parse_header(header):
    """
    Parse FASTA header to extract metadata.
    
    Format: GCF_002897375.1|WP_103129364.1|CDS|NZ_BFAG01000006.1:207150-207878(+)
    Returns: (full_header, genome_accession, protein_id, feature_type, seqid, start, end, strand)
    """
    parts = header.split('|')
    if len(parts) < 4:
        return None
    
    genome_acc = parts[0]
    protein_id = parts[1]
    feature_type = parts[2]
    coords_str = parts[3]
    
    # Parse coordinates: NZ_BFAG01000006.1:207150-207878(+)
    match = re.match(r'([^:]+):(\d+)-(\d+)\((.)\)', coords_str)
    if not match:
        return None
    
    seqid, start, end, strand = match.groups()
    
    return {
        'full_header': header,
        'protein_id': protein_id,
        'genome_accession': genome_acc,
        'feature_type': feature_type,
        'seqid': seqid,
        'start': int(start),
        'end': int(end),
        'strand': strand
    }

def parse_fasta(fasta_file):
    """
    Parse FASTA file, return list of (sequence_record_dict, sequence) tuples.
    Extracts metadata from headers.
    """
    sequences = []
    current_header = None
    current_seq = []
    metadata_rows = []
    
    with open(fasta_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header is not None:
                    header_data = parse_header(current_header)
                    if header_data:
                        sequences.append((header_data, ''.join(current_seq)))
                        metadata_rows.append(header_data)
                current_header = line[1:]  # Remove '>'
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_header is not None:
            header_data = parse_header(current_header)
            if header_data:
                sequences.append((header_data, ''.join(current_seq)))
                metadata_rows.append(header_data)
    
    metadata_df = pd.DataFrame(metadata_rows)
    return sequences, metadata_df

def get_embeddings(sequences, model, alphabet, device, batch_size=32):
    """
    Generate ESM2 embeddings for sequences.
    
    Args:
        sequences: list of (header_dict, sequence) tuples
        model: ESM2 model
        alphabet: alphabet for encoding
        device: torch device (cuda or cpu)
        batch_size: batch size for inference
    
    Returns:
        dict mapping protein_id → embedding (numpy array)
    """
    embeddings = {}
    
    # Process in batches
    for batch_start in tqdm(range(0, len(sequences), batch_size), desc="Embedding proteins"):
        batch_seqs = sequences[batch_start:batch_start + batch_size]
        
        # Prepare batch
        batch_headers = [seq_rec for seq_rec, _ in batch_seqs]
        batch_strs = [seq for _, seq in batch_seqs]
        batch_protein_ids = [h['protein_id'] for h in batch_headers]
        
        # Encode sequences
        batch_lens = [len(seq) for seq in batch_strs]
        
        # Create token batch (pad to max length in batch)
        max_len = max(batch_lens) + 2  # +2 for special tokens
        tokens = torch.zeros((len(batch_seqs), max_len), dtype=torch.int64, device=device)
        
        for i, seq in enumerate(batch_strs):
            seq_tokens = alphabet.encode(seq)
            tokens[i, :len(seq_tokens)] = torch.tensor(seq_tokens, device=device)
        
        # Get embeddings
        with torch.no_grad():
            results = model(tokens, repr_layers=[33])  # Layer 33 for ESM2-650M
        
        # Extract per-token embeddings, average over sequence length
        for i, protein_id in enumerate(batch_protein_ids):
            token_embeddings = results["representations"][33][i]  # (seq_len, embedding_dim)
            # Average over sequence (ignore padding)
            seq_len = batch_lens[i]
            avg_embedding = token_embeddings[:seq_len].mean(dim=0).cpu().numpy()
            embeddings[protein_id] = avg_embedding
    
    return embeddings

def main():
    parser = argparse.ArgumentParser(
        description="Generate ESM2 protein embeddings (per-family)"
    )
    parser.add_argument(
        "--fasta",
        required=True,
        help="Combined protein FASTA file for a single family"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output .npz file with embeddings"
    )
    parser.add_argument(
        '--output-metadata',
        required=True,
        help='path to output metadata file'
    )
    parser.add_argument(
        "--device",
        default="cuda" if torch.cuda.is_available() else "cpu",
        help="torch device (cuda or cpu)"
    )
    parser.add_argument(
        "--batch_size",
        type=int,
        default=32,
        help="Batch size for inference"
    )
    parser.add_argument(
        "--model",
        default="esm2_t33_650M_UR50D",
        help="ESM2 model name"
    )
    
    args = parser.parse_args()
    
    print(f"Loading ESM2 model: {args.model}")
    model, alphabet = esm.pretrained.load_model_and_alphabet_local(args.model)
    model = model.to(args.device)
    model.eval()
    
    print(f"Parsing FASTA: {args.fasta}")
    sequences, metadata_df = parse_fasta(args.fasta)
    print(f"Found {len(sequences)} sequences")
    print(f"Genomes: {metadata_df['genome_accession'].nunique()}")
    print(f"Unique genomes: {metadata_df['genome_accession'].unique()}")
    
    print(f"\nGenerating embeddings on {args.device}...")
    embeddings = get_embeddings(sequences, model, alphabet, args.device, args.batch_size)
    
    # Convert to numpy array (maintain protein ID order)
    protein_ids = list(embeddings.keys())
    embedding_matrix = np.array([embeddings[pid] for pid in protein_ids])
    
    print(f"Embedding shape: {embedding_matrix.shape}")
    
    # Create output directory
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Save embeddings, protein IDs, and metadata
    print(f"Saving to {args.output}")
    np.savez(args.output, 
             embeddings=embedding_matrix,
             protein_ids=np.array(protein_ids))
    
    # Also save metadata TSV alongside npz
    metadata_output = Path(args.output_metadata)
    metadata_df.to_csv(metadata_output, sep='\t', index=False)
    print(f"Saved metadata to {metadata_output}")
    
    print("Done!")

if __name__ == "__main__":
    main()
