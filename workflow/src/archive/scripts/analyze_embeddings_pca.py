#!/usr/bin/env python3
"""
PCA analysis of ESM2 protein embeddings (per-family) colored by genome accession.

Usage:
    python analyze_embeddings_pca.py \
        --embeddings results/00_controls/rfam/RF00050/esm_analysis/embeddings.npz \
        --output results/00_controls/rfam/RF00050/esm_analysis/pca_analysis/
"""

import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns

def load_embeddings_and_metadata(npz_file):
    """Load embeddings and protein IDs from npz file, load metadata from TSV."""
    # Load embeddings
    data = np.load(npz_file)
    embeddings = data['embeddings']
    protein_ids = data['protein_ids']
    
    # Load metadata TSV from same directory
    metadata_file = Path(npz_file).parent / "metadata.tsv"
    if metadata_file.exists():
        metadata = pd.read_csv(metadata_file, sep='\t')
    else:
        print(f"Warning: metadata.tsv not found at {metadata_file}")
        # Create minimal metadata from protein IDs
        metadata = pd.DataFrame({'protein_id': protein_ids})
    
    return protein_ids, embeddings, metadata

def compute_pca(embeddings, n_components=50):
    """Compute PCA on embeddings."""
    # Standardize
    scaler = StandardScaler()
    embeddings_scaled = scaler.fit_transform(embeddings)
    
    # PCA
    pca = PCA(n_components=n_components)
    pca_coords = pca.fit_transform(embeddings_scaled)
    
    print(f"PCA explained variance ratio (first 10 components):")
    for i, var in enumerate(pca.explained_variance_ratio_[:10]):
        print(f"  PC{i+1}: {var:.2%}")
    
    return pca_coords, pca

def plot_pca_2d(pca_coords, protein_ids, metadata, color_by, output_dir):
    """Plot 2D PCA colored by metadata."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create DataFrame with PCA coords
    df_plot = pd.DataFrame({
        'protein_id': protein_ids,
        'PC1': pca_coords[:, 0],
        'PC2': pca_coords[:, 1]
    })
    
    # Merge with metadata
    df_plot = df_plot.merge(metadata, on='protein_id', how='left')
    
    # Plot
    fig, ax = plt.subplots(figsize=(12, 9))
    
    # Get unique values and assign colors
    unique_vals = df_plot[color_by].dropna().unique()
    n_colors = len(unique_vals)
    
    # Use palette if reasonable number of categories
    if n_colors <= 20:
        palette = sns.color_palette("husl", n_colors)
        sns.scatterplot(data=df_plot, x='PC1', y='PC2', 
                       hue=color_by, palette=palette, 
                       s=100, alpha=0.7, ax=ax)
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8, title=color_by)
    else:
        # Too many categories, use index
        color_map = {val: i for i, val in enumerate(unique_vals)}
        colors = df_plot[color_by].map(color_map)
        scatter = ax.scatter(df_plot['PC1'], df_plot['PC2'], 
                            c=colors, s=100, alpha=0.7, cmap='tab20')
        plt.colorbar(scatter, ax=ax, label=color_by)
    
    ax.set_xlabel(f'PC1')
    ax.set_ylabel(f'PC2')
    ax.set_title(f'ESM2 Protein Embeddings PCA (colored by {color_by})')
    ax.grid(True, alpha=0.3)
    
    output_file = output_dir / f"pca_2d_{color_by}.png"
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()
    
    return df_plot

def plot_pca_3d(pca_coords, protein_ids, metadata, color_by, output_dir):
    """Plot 3D PCA colored by metadata."""
    if pca_coords.shape[1] < 3:
        print("Warning: Need at least 3 components for 3D plot, skipping")
        return
    
    output_dir = Path(output_dir)
    
    from mpl_toolkits.mplot3d import Axes3D
    
    df_plot = pd.DataFrame({
        'protein_id': protein_ids,
        'PC1': pca_coords[:, 0],
        'PC2': pca_coords[:, 1],
        'PC3': pca_coords[:, 2]
    })
    
    df_plot = df_plot.merge(metadata, on='protein_id', how='left')
    
    fig = plt.figure(figsize=(14, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    unique_vals = df_plot[color_by].dropna().unique()
    palette = sns.color_palette("husl", len(unique_vals))
    
    for val in unique_vals:
        mask = df_plot[color_by] == val
        ax.scatter(df_plot[mask]['PC1'], df_plot[mask]['PC2'], df_plot[mask]['PC3'],
                  label=str(val)[:30], s=50, alpha=0.7)  # Truncate long labels
    
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    ax.set_zlabel('PC3')
    ax.set_title(f'ESM2 Protein Embeddings 3D PCA (colored by {color_by})')
    ax.legend(fontsize=8)
    
    output_file = output_dir / f"pca_3d_{color_by}.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()

def plot_variance_explained(pca, output_dir, n_components=20):
    """Plot cumulative variance explained."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    cumvar = np.cumsum(pca.explained_variance_ratio_[:n_components])
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(range(1, n_components+1), cumvar, 'bo-', linewidth=2, markersize=6)
    ax.axhline(y=0.95, color='r', linestyle='--', label='95% variance')
    ax.set_xlabel('Number of Components')
    ax.set_ylabel('Cumulative Explained Variance')
    ax.set_title('PCA Variance Explained')
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    output_file = output_dir / "variance_explained.png"
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    print(f"Saved: {output_file}")
    plt.close()

def write_pca_table(df_plot, output_dir):
    """Write PCA coordinates + metadata to TSV."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / "pca_coordinates.tsv"
    df_plot.to_csv(output_file, sep='\t', index=False)
    print(f"Saved: {output_file}")

def main():
    parser = argparse.ArgumentParser(
        description="Analyze ESM2 embeddings with PCA (per-family)"
    )
    parser.add_argument("--embeddings", required=True, help="embeddings.npz file")
    parser.add_argument("--color_by", default="genome_accession",
                       help="Metadata column to color by")
    parser.add_argument("--output", required=True,
                       help="Output directory")
    parser.add_argument("--n_components_pca", type=int, default=50,
                       help="Number of PCA components to compute")
    
    args = parser.parse_args()
    
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load data
    print(f"Loading embeddings from {args.embeddings}")
    protein_ids, embeddings, metadata = load_embeddings_and_metadata(args.embeddings)
    
    print(f"Loaded {len(protein_ids)} proteins")
    print(f"Metadata columns: {list(metadata.columns)}")
    
    if args.color_by not in metadata.columns:
        print(f"Error: '{args.color_by}' not found in metadata columns")
        print(f"Available columns: {list(metadata.columns)}")
        exit(1)
    
    print(f"Unique {args.color_by}: {metadata[args.color_by].nunique()}")
    print(f"  Values: {metadata[args.color_by].unique()}")
    
    # Compute PCA
    print(f"\nComputing PCA with {args.n_components_pca} components...")
    pca_coords, pca = compute_pca(embeddings, n_components=args.n_components_pca)
    
    # Visualizations
    print("\nGenerating plots...")
    
    # 2D PCA
    df_plot = plot_pca_2d(pca_coords, protein_ids, metadata, 
                          args.color_by, output_dir)
    
    # 3D PCA (if enough components)
    if args.n_components_pca >= 3:
        plot_pca_3d(pca_coords, protein_ids, metadata, 
                   args.color_by, output_dir)
    
    # Variance explained
    plot_variance_explained(pca, output_dir, n_components=min(20, args.n_components_pca))
    
    # Write PCA table
    write_pca_table(df_plot, output_dir)
    
    print(f"\nAll outputs written to {output_dir}")

if __name__ == "__main__":
    main()
