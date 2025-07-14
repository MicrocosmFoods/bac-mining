#!/usr/bin/env python3
"""
Create comprehensive genome summaries by combining results from bac-mining workflow.

This script combines the main results from the bac-mining workflow to create
per-genome summaries with locus tag information, annotations, and peptide/BGC types.

Usage:
    python create_genome_summaries.py --input-dir /path/to/main_results --output-dir /path/to/output
"""

import argparse
import os
import sys
import pandas as pd
from pathlib import Path
import glob


def find_result_files(input_dir):
    """Find all the result files in the main_results directory."""
    input_path = Path(input_dir)
    
    # Find the main result files
    files = {
        'smorfinder': input_path / 'all_smorfinder_results.tsv',
        'deeppeptide': input_path / 'all_deeppeptide_results.tsv', 
        'antismash': input_path / 'bgc_info' / 'antismash_summary.tsv',
        'kofamscan': input_path / 'combined_kofamscan_results.tsv'
    }
    
    # Check which files exist
    existing_files = {}
    for name, filepath in files.items():
        if filepath.exists():
            existing_files[name] = filepath
            print(f"Found {name} results: {filepath}")
        else:
            print(f"Warning: {name} results not found at {filepath}")
    
    return existing_files


def create_genome_summaries(input_dir, output_dir, genome_list=None):
    """Create genome summaries by combining all available results."""
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Find result files
    result_files = find_result_files(input_dir)
    
    if not result_files:
        print("Error: No result files found!")
        sys.exit(1)
    
    # Load and process each result file
    results = {}
    
    # Load smorfinder results if available
    if 'smorfinder' in result_files:
        print("Loading smorfinder results...")
        smorf_df = pd.read_csv(result_files['smorfinder'], sep='\t')
        if not smorf_df.empty:
            results['smorfinder'] = smorf_df
            print(f"Loaded {len(smorf_df)} smorfinder entries")
    
    # Load deeppeptide results if available
    if 'deeppeptide' in result_files:
        print("Loading deeppeptide results...")
        deepp_df = pd.read_csv(result_files['deeppeptide'], sep='\t')
        if not deepp_df.empty:
            results['deeppeptide'] = deepp_df
            print(f"Loaded {len(deepp_df)} deeppeptide entries")
    
    # Load antismash results if available
    if 'antismash' in result_files:
        print("Loading antismash results...")
        antismash_df = pd.read_csv(result_files['antismash'], sep='\t')
        if not antismash_df.empty:
            results['antismash'] = antismash_df
            print(f"Loaded {len(antismash_df)} antismash entries")
    
    # Load kofamscan results if available
    if 'kofamscan' in result_files:
        print("Loading kofamscan results...")
        kofam_df = pd.read_csv(result_files['kofamscan'], sep='\t')
        if not kofam_df.empty:
            results['kofamscan'] = kofam_df
            print(f"Loaded {len(kofam_df)} kofamscan entries")
    
    # Get list of genomes to process
    if genome_list:
        with open(genome_list, 'r') as f:
            genomes = [line.strip() for line in f if line.strip()]
    else:
        # Extract genome names from any available result file
        genomes = set()
        for name, df in results.items():
            if 'genome' in df.columns:
                genomes.update(df['genome'].unique())
        genomes = sorted(list(genomes))
    
    print(f"Processing {len(genomes)} genomes...")
    
    # Process each genome
    for genome in genomes:
        print(f"Processing genome: {genome}")
        
        # Initialize summary data for this genome
        summary_data = []
        
        # Add smorfinder results
        if 'smorfinder' in results:
            smorf_genome = results['smorfinder'][results['smorfinder']['genome'] == genome]
            for _, row in smorf_genome.iterrows():
                summary_data.append({
                    'genome': genome,
                    'locus_tag': row.get('locus_tag', ''),
                    'start': row.get('start', ''),
                    'stop': row.get('stop', ''),
                    'peptide_type': 'smorfinder',
                    'peptide_name': row.get('peptide_name', ''),
                    'bgc_type': '',
                    'kegg_annotation': '',
                    'source_file': 'smorfinder'
                })
        
        # Add deeppeptide results
        if 'deeppeptide' in results:
            deepp_genome = results['deeppeptide'][results['deeppeptide']['genome'] == genome]
            for _, row in deepp_genome.iterrows():
                summary_data.append({
                    'genome': genome,
                    'locus_tag': row.get('locus_tag', ''),
                    'start': row.get('start', ''),
                    'stop': row.get('stop', ''),
                    'peptide_type': 'deeppeptide',
                    'peptide_name': row.get('peptide_name', ''),
                    'bgc_type': '',
                    'kegg_annotation': '',
                    'source_file': 'deeppeptide'
                })
        
        # Add antismash results
        if 'antismash' in results:
            antismash_genome = results['antismash'][results['antismash']['genome'] == genome]
            for _, row in antismash_genome.iterrows():
                summary_data.append({
                    'genome': genome,
                    'locus_tag': row.get('locus_tag', ''),
                    'start': row.get('start', ''),
                    'stop': row.get('stop', ''),
                    'peptide_type': '',
                    'peptide_name': '',
                    'bgc_type': row.get('bgc_type', ''),
                    'kegg_annotation': '',
                    'source_file': 'antismash'
                })
        
        # Add kofamscan results
        if 'kofamscan' in results:
            kofam_genome = results['kofamscan'][results['kofamscan']['genome'] == genome]
            for _, row in kofam_genome.iterrows():
                summary_data.append({
                    'genome': genome,
                    'locus_tag': row.get('locus_tag', ''),
                    'start': '',
                    'stop': '',
                    'peptide_type': '',
                    'peptide_name': '',
                    'bgc_type': '',
                    'kegg_annotation': row.get('annotation', ''),
                    'source_file': 'kofamscan'
                })
        
        # Create summary DataFrame for this genome
        if summary_data:
            summary_df = pd.DataFrame(summary_data)
            
            # Save genome summary
            output_file = output_path / f"{genome}_genome_summary.tsv"
            summary_df.to_csv(output_file, sep='\t', index=False)
            print(f"Saved {len(summary_df)} entries to {output_file}")
        else:
            print(f"No data found for genome {genome}")
    
    # Create combined summary
    all_summaries = []
    for genome in genomes:
        genome_file = output_path / f"{genome}_genome_summary.tsv"
        if genome_file.exists():
            genome_df = pd.read_csv(genome_file, sep='\t')
            all_summaries.append(genome_df)
    
    if all_summaries:
        combined_df = pd.concat(all_summaries, ignore_index=True)
        combined_file = output_path / "all_genome_summaries.tsv"
        combined_df.to_csv(combined_file, sep='\t', index=False)
        print(f"Saved combined summary with {len(combined_df)} entries to {combined_file}")
    
    print("Genome summary creation complete!")


def main():
    parser = argparse.ArgumentParser(
        description="Create comprehensive genome summaries from bac-mining results"
    )
    parser.add_argument(
        "--input-dir", 
        required=True,
        help="Path to main_results directory from bac-mining workflow"
    )
    parser.add_argument(
        "--output-dir", 
        required=True,
        help="Output directory for genome summaries"
    )
    parser.add_argument(
        "--genome-list",
        help="Optional file containing list of genomes to process (one per line)"
    )
    
    args = parser.parse_args()
    
    # Validate input directory
    if not os.path.exists(args.input_dir):
        print(f"Error: Input directory {args.input_dir} does not exist!")
        sys.exit(1)
    
    # Create summaries
    create_genome_summaries(args.input_dir, args.output_dir, args.genome_list)


if __name__ == "__main__":
    main() 