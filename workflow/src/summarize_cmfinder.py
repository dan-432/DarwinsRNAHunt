
"""
Summarize CMfinder analysis results
"""
import argparse
from pathlib import Path
import re
from datetime import datetime

def read_file_safe(filepath):
    """Safely read a file, return empty string if not found."""
    try:
        with open(filepath, 'r') as f:
            return f.read()
    except FileNotFoundError:
        return "File not found"

def count_hits(tblout_file):
    """Count number of hits in a cmsearch tblout file."""
    try:
        with open(tblout_file, 'r') as f:
            # Count non-comment lines
            return sum(1 for line in f if line.strip() and not line.startswith('#'))
    except FileNotFoundError:
        return 0

def get_motif_stats(cmfinder_dir):
    """Collect statistics for all motifs."""
    cmfinder_path = Path(cmfinder_dir)
    
    # Find all final models
    final_models_dir = cmfinder_path / "02_models" / "final"
    if not final_models_dir.exists():
        print(f"WARNING: Final models directory not found: {final_models_dir}")
        return []
    
    final_models = list(final_models_dir.glob("*_final.cm"))
    
    print(f"Found {len(final_models)} final models in {final_models_dir}")
    
    motif_stats = []
    for model in final_models:
        # Extract motif ID from filename like "h1_1_final.cm"
        motif_id = model.stem.replace('_final', '')
        
        print(f"Processing motif: {motif_id}")
        
        # Count search hits
        hits_file = cmfinder_path / "03_homolog_search" / "hits" / f"{motif_id}_hits.tblout"
        hit_count = count_hits(hits_file)
        
        print(f"  Hits file: {hits_file}, count: {hit_count}")
        
        # Check if visualization exists
        viz_file = cmfinder_path / "05_visualizations" / f"{motif_id}_structure.pdf"
        has_viz = viz_file.exists() and viz_file.stat().st_size > 0
        
        # Get model info from build log
        build_log = cmfinder_path.parent.parent / "logs" / "cmfinder" / f"*_02_build_cm_{motif_id}.log"
        
        motif_stats.append({
            'id': motif_id,
            'hits': hit_count,
            'has_visualization': has_viz,
            'model_path': str(model)
        })
    
    return motif_stats

def summarize_cmfinder(target_domain, cmfinder_dir, output_summary, output_report):
    """Generate summary of CMfinder analysis."""
    
    cmfinder_path = Path(cmfinder_dir)
    
    print(f"Summarizing CMfinder results for {target_domain}")
    print(f"CMfinder directory: {cmfinder_path}")
    
    # Get motif statistics
    motif_stats = get_motif_stats(cmfinder_dir)
    
    if not motif_stats:
        print("WARNING: No motif statistics found!")
    
    # Collect overall statistics
    stats = {
        'target_domain': target_domain,
        'timestamp': datetime.now().isoformat(),
        'total_motifs': len(motif_stats),
        'total_hits': sum(m['hits'] for m in motif_stats),
        'with_visualizations': sum(1 for m in motif_stats if m['has_visualization'])
    }
    
    print(f"Overall stats: {stats}")
    
    # Write text summary
    with open(output_summary, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write(f"CMfinder Analysis Summary: {target_domain}\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Analysis Date: {stats['timestamp']}\n")
        f.write(f"Pipeline: Based on Narunsky et al., NAR 2024\n")
        f.write(f"DOI: 10.1093/nar/gkae248\n\n")
        
        f.write("Overall Statistics:\n")
        f.write("-" * 40 + "\n")
        f.write(f"  Total motifs discovered:    {stats['total_motifs']}\n")
        f.write(f"  Total homolog hits found:   {stats['total_hits']}\n")
        f.write(f"  Motifs with visualization:  {stats['with_visualizations']}\n\n")
        
        if motif_stats:
            f.write("Individual Motif Results:\n")
            f.write("-" * 40 + "\n")
            for motif in sorted(motif_stats, key=lambda x: x['id']):
                f.write(f"\n  Motif: {motif['id']}\n")
                f.write(f"    Homolog hits: {motif['hits']}\n")
                f.write(f"    Visualization: {'Yes' if motif['has_visualization'] else 'No'}\n")
                f.write(f"    Model: {motif['model_path']}\n")
        else:
            f.write("No motifs found!\n")
        
        f.write("\n" + "=" * 40 + "\n")
        f.write("Output Directories:\n")
        f.write("-" * 40 + "\n")
        f.write(f"  Initial results:     {cmfinder_path / '01_initial_discovery'}\n")
        f.write(f"  Final models:        {cmfinder_path / '02_models' / 'final'}\n")
        f.write(f"  Visualizations:      {cmfinder_path / '05_visualizations'}\n")
        f.write(f"  Search results:      {cmfinder_path / '03_homolog_search' / 'hits'}\n")
    
    # Write HTML report
    if motif_stats:
        motif_rows = "\n".join([
            f"<tr><td>{m['id']}</td><td>{m['hits']}</td><td>{'✓' if m['has_visualization'] else '✗'}</td></tr>"
            for m in sorted(motif_stats, key=lambda x: x['id'])
        ])
    else:
        motif_rows = "<tr><td colspan='3'>No motifs found</td></tr>"
    
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>CMfinder Analysis Report: {target_domain}</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 40px; background-color: #f5f5f5; }}
            .container {{ max-width: 1200px; margin: 0 auto; background-color: white; padding: 30px; box-shadow: 0 0 10px rgba(0,0,0,0.1); }}
            h1 {{ color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }}
            h2 {{ color: #34495e; border-bottom: 2px solid #3498db; padding-bottom: 5px; margin-top: 30px; }}
            table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
            th, td {{ border: 1px solid #ddd; padding: 12px; text-align: left; }}
            th {{ background-color: #3498db; color: white; }}
            tr:nth-child(even) {{ background-color: #f2f2f2; }}
            tr:hover {{ background-color: #e8f4f8; }}
            .summary-box {{ background-color: #ecf0f1; padding: 20px; border-radius: 5px; margin: 20px 0; border-left: 4px solid #3498db; }}
            .stats {{ display: grid; grid-template-columns: repeat(3, 1fr); gap: 20px; margin: 20px 0; }}
            .stat-card {{ background-color: #3498db; color: white; padding: 20px; border-radius: 5px; text-align: center; }}
            .stat-number {{ font-size: 2.5em; font-weight: bold; }}
            .stat-label {{ font-size: 1em; margin-top: 5px; }}
            .warning {{ background-color: #fff3cd; border: 1px solid #ffc107; padding: 15px; margin: 20px 0; border-radius: 5px; }}
        </style>
    </head>
    <body>
        <div class="container">
            <h1>CMfinder RNA Structural Motif Analysis</h1>
            <h2>Target Domain: {target_domain}</h2>
            
            <div class="summary-box">
                <strong>Analysis Date:</strong> {stats['timestamp']}<br>
                <strong>Pipeline:</strong> Based on Narunsky et al., Nucleic Acids Research, 2024<br>
                <strong>DOI:</strong> 10.1093/nar/gkae248
            </div>
            
            <h2>Overall Statistics</h2>
            <div class="stats">
                <div class="stat-card">
                    <div class="stat-number">{stats['total_motifs']}</div>
                    <div class="stat-label">Motifs Discovered</div>
                </div>
                <div class="stat-card">
                    <div class="stat-number">{stats['total_hits']}</div>
                    <div class="stat-label">Total Homolog Hits</div>
                </div>
                <div class="stat-card">
                    <div class="stat-number">{stats['with_visualizations']}</div>
                    <div class="stat-label">Visualizations</div>
                </div>
            </div>
            
            <h2>Individual Motif Results</h2>
            <table>
                <tr>
                    <th>Motif ID</th>
                    <th>Homolog Hits</th>
                    <th>Visualization</th>
                </tr>
                {motif_rows}
            </table>
            
            <h2>Motif Details</h2>
            <div class="summary-box">
                <strong>h1 motifs:</strong> Single stem-loop structures (simpler RNA structures)<br>
                <strong>h2 motifs:</strong> Double stem-loop structures (more complex RNA structures)<br>
                <strong>Note:</strong> Numbers like h1_1, h1_2 indicate the rank (best to worst by CMfinder score)
            </div>
            
            <h2>Output Directories</h2>
            <ul>
                <li><strong>Initial CMfinder results:</strong> {cmfinder_path / '01_initial_discovery'}</li>
                <li><strong>Built models:</strong> {cmfinder_path / '02_models' / 'built'}</li>
                <li><strong>Calibrated models:</strong> {cmfinder_path / '02_models' / 'calibrated'}</li>
                <li><strong>Final models:</strong> {cmfinder_path / '02_models' / 'final'}</li>
                <li><strong>Search results:</strong> {cmfinder_path / '03_homolog_search' / 'hits'}</li>
                <li><strong>Expanded sequences:</strong> {cmfinder_path / '03_homolog_search' / 'expanded_seqs'}</li>
                <li><strong>Refined alignments:</strong> {cmfinder_path / '04_refinement'}</li>
                <li><strong>Structure visualizations:</strong> {cmfinder_path / '05_visualizations'}</li>
            </ul>
            
            <h2>Interpretation</h2>
            <div class="summary-box">
                <p><strong>Your Results:</strong> CMfinder discovered {stats['total_motifs']} high-quality structural RNA motifs 
                with {stats['total_hits']} total homolog hits across your sequences.</p>
                
                <p><strong>What this means:</strong></p>
                <ul>
                    <li>Each motif represents a conserved RNA structure found in your intergenic regions</li>
                    <li>These could be regulatory elements like riboswitches, RNA thermometers, or other functional RNAs</li>
                    <li>Hits indicate how many of your input sequences contain each motif</li>
                    <li>Higher hit counts suggest more widespread conservation</li>
                </ul>
            </div>
        </div>
    </body>
    </html>
    """
    
    with open(output_report, 'w') as f:
        f.write(html_content)
    
    print(f"Summary written to: {output_summary}")
    print(f"HTML report written to: {output_report}")
    print(f"\nFound {stats['total_motifs']} motifs with {stats['total_hits']} total homolog hits")

def main():
    parser = argparse.ArgumentParser(
        description='Summarize CMfinder analysis results'
    )
    parser.add_argument('--target-domain', required=True)
    parser.add_argument('--cmfinder-dir', required=True)
    parser.add_argument('--output-summary', required=True)
    parser.add_argument('--output-report', required=True)
    
    args = parser.parse_args()
    summarize_cmfinder(
        args.target_domain,
        args.cmfinder_dir,
        args.output_summary,
        args.output_report
    )

if __name__ == '__main__':
    main()