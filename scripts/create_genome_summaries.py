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


def find_ffn_files(predicted_orfs_dir):
    """Find all .ffn files in the predicted_orfs directory."""
    ffn_path = Path(predicted_orfs_dir)
    
    if not ffn_path.exists():
        print(f"Warning: Predicted ORFs directory {predicted_orfs_dir} not found!")
        return {}
    
    # Find all .ffn files
    ffn_files = {}
    for ffn_file in ffn_path.glob("*.ffn"):
        genome_name = ffn_file.stem  # Remove .ffn extension
        ffn_files[genome_name] = ffn_file
        print(f"Found .ffn file for genome {genome_name}: {ffn_file}")
    
    return ffn_files


def parse_ffn_file(ffn_file):
    """Parse .ffn file to extract gene information."""
    genes = {}
    
    with open(ffn_file, 'r') as f:
        current_gene = None
        current_sequence = ""
        
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous gene if exists
                if current_gene:
                    genes[current_gene] = {
                        'sequence': current_sequence,
                        'length': len(current_sequence)
                    }
                
                # Parse header line
                header = line[1:]  # Remove '>'
                parts = header.split()
                current_gene = parts[0]  # First part is usually the locus tag
                current_sequence = ""
            else:
                current_sequence += line
        
        # Don't forget the last gene
        if current_gene:
            genes[current_gene] = {
                'sequence': current_sequence,
                'length': len(current_sequence)
            }
    
    return genes


def parse_ffn_file_with_coordinates(ffn_file):
    """Parse .ffn file to extract gene information with coordinates from header."""
    genes = {}
    
    with open(ffn_file, 'r') as f:
        current_gene = None
        current_sequence = ""
        current_start = None
        current_end = None
        
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous gene if exists
                if current_gene:
                    genes[current_gene] = {
                        'sequence': current_sequence,
                        'length': len(current_sequence),
                        'start': current_start,
                        'end': current_end
                    }
                
                # Parse header line
                header = line[1:]  # Remove '>'
                parts = header.split()
                current_gene = parts[0]  # First part is the locus tag (clean, no whitespace)
                current_sequence = ""
                current_start = None
                current_end = None
                
                # Try to extract coordinates from header description
                if len(parts) > 1:
                    description = ' '.join(parts[1:])
                    # Look for [location=start..end] pattern
                    if '[location=' in description:
                        location_part = description.split('[location=')[1].split(']')[0]
                        if '..' in location_part:
                            try:
                                start_str, end_str = location_part.split('..')
                                current_start = int(start_str)
                                current_end = int(end_str)
                            except ValueError:
                                pass
            else:
                current_sequence += line
        
        # Don't forget the last gene
        if current_gene:
            genes[current_gene] = {
                'sequence': current_sequence,
                'length': len(current_sequence),
                'start': current_start,
                'end': current_end
            }
    
    return genes


def get_genome_column_name(df, file_type):
    """Determine the correct column name for genome identifier in a DataFrame."""
    possible_names = ['genome', 'genome_name', 'mag_id', 'sample', 'strain']
    
    for col_name in possible_names:
        if col_name in df.columns:
            return col_name
    
    # If none found, print available columns and return None
    print(f"Warning: Could not find genome column in {file_type} file. Available columns: {list(df.columns)}")
    return None


def get_locus_column_name(df, file_type):
    """Determine the correct column name for locus tag in a DataFrame."""
    possible_names = ['locus_tag', 'gene_name', 'gene_id', 'orf_id', 'id']
    
    for col_name in possible_names:
        if col_name in df.columns:
            return col_name
    
    # If none found, print available columns and return None
    print(f"Warning: Could not find locus tag column in {file_type} file. Available columns: {list(df.columns)}")
    return None


def create_genome_summaries(input_dir, output_dir, predicted_orfs_dir, genome_list=None):
    """Create genome summaries by combining all available results with .ffn gene information."""
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Find result files
    result_files = find_result_files(input_dir)
    
    # Find .ffn files
    ffn_files = find_ffn_files(predicted_orfs_dir)
    
    if not result_files and not ffn_files:
        print("Error: No result files or .ffn files found!")
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
            print(f"  Columns: {list(smorf_df.columns)}")
    
    # Load deeppeptide results if available
    if 'deeppeptide' in result_files:
        print("Loading deeppeptide results...")
        deepp_df = pd.read_csv(result_files['deeppeptide'], sep='\t')
        if not deepp_df.empty:
            results['deeppeptide'] = deepp_df
            print(f"Loaded {len(deepp_df)} deeppeptide entries")
            print(f"  Columns: {list(deepp_df.columns)}")
    
    # Load antismash results if available
    if 'antismash' in result_files:
        print("Loading antismash results...")
        antismash_df = pd.read_csv(result_files['antismash'], sep='\t')
        if not antismash_df.empty:
            results['antismash'] = antismash_df
            print(f"Loaded {len(antismash_df)} antismash entries")
            print(f"  Columns: {list(antismash_df.columns)}")
    
    # Load kofamscan results if available
    if 'kofamscan' in result_files:
        print("Loading kofamscan results...")
        kofam_df = pd.read_csv(result_files['kofamscan'], sep='\t')
        if not kofam_df.empty:
            results['kofamscan'] = kofam_df
            print(f"Loaded {len(kofam_df)} kofamscan entries")
            print(f"  Columns: {list(kofam_df.columns)}")
            # Show first few rows to understand the data structure
            print(f"  First few rows:")
            print(kofam_df.head(3).to_string())
    
    # Get list of genomes to process
    if genome_list:
        with open(genome_list, 'r') as f:
            genomes = [line.strip() for line in f if line.strip()]
    else:
        # Use genomes from .ffn files, or fall back to result files
        if ffn_files:
            genomes = sorted(list(ffn_files.keys()))
        else:
            # Extract genome names from any available result file
            genomes = set()
            for name, df in results.items():
                genome_col = get_genome_column_name(df, name)
                if genome_col:
                    genomes.update(df[genome_col].unique())
            genomes = sorted(list(genomes))
    
    print(f"Processing {len(genomes)} genomes...")
    
    # Process each genome
    all_summaries = []
    
    for genome in genomes:
        print(f"Processing genome: {genome}")
        
        # Parse .ffn file if available
        genes = {}
        if genome in ffn_files:
            print(f"  Parsing .ffn file for {genome}...")
            genes = parse_ffn_file_with_coordinates(ffn_files[genome])
            print(f"  Found {len(genes)} genes in .ffn file")
        
        # Initialize summary data for this genome
        summary_data = []
        
        # Create entries for all genes from .ffn file
        for locus_tag, gene_info in genes.items():
            summary_data.append({
                'genome': genome,
                'locus_tag': locus_tag,
                'gene_length': gene_info['length'],
                'smorfinder_peptide': '',
                'smorfinder_peptide_name': '',
                'deeppeptide_peptide': '',
                'deeppeptide_peptide_name': '',
                'antismash_bgc_type': '',
                'ko_number': '',
                'kofamscan_annotation': ''
            })
        
        # Add annotations from result files
        annotations_added = 0
        
        # Add smorfinder results
        if 'smorfinder' in results:
            smorf_df = results['smorfinder']
            genome_col = get_genome_column_name(smorf_df, 'smorfinder')
            
            if genome_col:
                smorf_genome = smorf_df[smorf_df[genome_col] == genome]
                for _, row in smorf_genome.iterrows():
                    locus_tag = row.get('seqid', '')  # Use seqid as locus tag
                    smorf_start = row.get('start', None)
                    smorf_end = row.get('end', None)
                    
                    # Try to match by locus tag first
                    matched = False
                    for entry in summary_data:
                        if entry['locus_tag'] == locus_tag:
                            entry['smorfinder_peptide'] = 'smorfinder'
                            entry['smorfinder_peptide_name'] = row.get('smorfam', '')  # Use smorfam as peptide name
                            annotations_added += 1
                            matched = True
                            break
                    
                    # If no locus tag match, try coordinate matching (if coordinates are available)
                    if not matched and smorf_start and smorf_end:
                        for entry in summary_data:
                            if (hasattr(entry, 'gene_start') and hasattr(entry, 'gene_end') and
                                entry['gene_start'] and entry['gene_end'] and
                                entry['gene_start'] == smorf_start and entry['gene_end'] == smorf_end):
                                entry['smorfinder_peptide'] = 'smorfinder'
                                entry['smorfinder_peptide_name'] = row.get('smorfam', '')
                                annotations_added += 1
                                matched = True
                                break
                    
                    if not matched:
                        # Gene not in .ffn file, add new entry
                        summary_data.append({
                            'genome': genome,
                            'locus_tag': locus_tag if locus_tag else f"smorf_{smorf_start}_{smorf_end}",
                            'gene_length': '',
                            'smorfinder_peptide': 'smorfinder',
                            'smorfinder_peptide_name': row.get('smorfam', ''),
                            'deeppeptide_peptide': '',
                            'deeppeptide_peptide_name': '',
                            'antismash_bgc_type': '',
                            'ko_number': '',
                            'kofamscan_annotation': ''
                        })
        
        # Add deeppeptide results
        if 'deeppeptide' in results:
            deepp_df = results['deeppeptide']
            genome_col = get_genome_column_name(deepp_df, 'deeppeptide')
            
            if genome_col:
                deepp_genome = deepp_df[deepp_df[genome_col] == genome]
                for _, row in deepp_genome.iterrows():
                    peptide_id = row.get('peptide_id', '')
                    
                    # Extract gene name from peptide_id (everything before _start)
                    if peptide_id and '_start' in peptide_id:
                        locus_tag = peptide_id.split('_start')[0]
                    else:
                        locus_tag = peptide_id
                    
                    # Try to match by locus tag
                    matched = False
                    for entry in summary_data:
                        if entry['locus_tag'] == locus_tag:
                            entry['deeppeptide_peptide'] = 'deeppeptide'
                            entry['deeppeptide_peptide_name'] = row.get('peptide_class', '')  # Use peptide_class as name
                            annotations_added += 1
                            matched = True
                            break
                    
                    if not matched:
                        # Gene not in .ffn file, add new entry
                        summary_data.append({
                            'genome': genome,
                            'locus_tag': locus_tag,
                            'gene_length': '',
                            'smorfinder_peptide': '',
                            'smorfinder_peptide_name': '',
                            'deeppeptide_peptide': 'deeppeptide',
                            'deeppeptide_peptide_name': row.get('peptide_class', ''),
                            'antismash_bgc_type': '',
                            'ko_number': '',
                            'kofamscan_annotation': ''
                        })
        
        # Add antismash results
        if 'antismash' in results:
            antismash_df = results['antismash']
            genome_col = get_genome_column_name(antismash_df, 'antismash')
            locus_col = get_locus_column_name(antismash_df, 'antismash')
            
            if genome_col:
                antismash_genome = antismash_df[antismash_df[genome_col] == genome]
                for _, row in antismash_genome.iterrows():
                    locus_tag = row.get(locus_col, '') if locus_col else ''
                    
                    if locus_tag:
                        # Try to match by locus tag first
                        matched = False
                        for entry in summary_data:
                            if entry['locus_tag'] == locus_tag:
                                entry['antismash_bgc_type'] = row.get('bgc_type', '')
                                annotations_added += 1
                                matched = True
                                break
                        
                        if not matched:
                            # Gene not in .ffn file, add new entry
                            summary_data.append({
                                'genome': genome,
                                'locus_tag': locus_tag,
                                'gene_length': '',
                                'smorfinder_peptide': '',
                                'smorfinder_peptide_name': '',
                                'deeppeptide_peptide': '',
                                'deeppeptide_peptide_name': '',
                                'antismash_bgc_type': row.get('bgc_type', ''),
                                'ko_number': '',
                                'kofamscan_annotation': ''
                            })
        
        # Add kofamscan results
        if 'kofamscan' in results:
            kofam_df = results['kofamscan']
            genome_col = get_genome_column_name(kofam_df, 'kofamscan')
            
            if genome_col:
                kofam_genome = kofam_df[kofam_df[genome_col] == genome]
                
                # Group kofamscan results by gene_name to find best scoring annotation per gene
                gene_annotations = {}
                for _, row in kofam_genome.iterrows():
                    locus_tag = row.get('gene_name', '')  # Use gene_name as locus tag
                    
                    if locus_tag:
                        ko_number = row.get('KO', '')  # KO number
                        definition = row.get('KO_definition', '')  # Annotation description
                        score = row.get('score', 0)  # Score for ranking
                        threshold = row.get('threshold', 0)  # Threshold score
                        
                        # Calculate confidence score (higher is better)
                        confidence_score = score - threshold if score and threshold else 0
                        
                        # Keep the annotation with highest confidence score
                        if locus_tag not in gene_annotations or confidence_score > gene_annotations[locus_tag]['confidence']:
                            gene_annotations[locus_tag] = {
                                'ko_number': ko_number,
                                'definition': definition,
                                'confidence': confidence_score
                            }
                
                # Now add the best annotations to summary data
                for locus_tag, annotation in gene_annotations.items():
                    matched = False
                    for entry in summary_data:
                        if entry['locus_tag'] == locus_tag:
                            entry['ko_number'] = annotation['ko_number']
                            entry['kofamscan_annotation'] = annotation['definition']
                            annotations_added += 1
                            matched = True
                            break
                    
                    if not matched:
                        # Gene not in .ffn file, add new entry
                        summary_data.append({
                            'genome': genome,
                            'locus_tag': locus_tag,
                            'gene_length': '',
                            'smorfinder_peptide': '',
                            'smorfinder_peptide_name': '',
                            'deeppeptide_peptide': '',
                            'deeppeptide_peptide_name': '',
                            'antismash_bgc_type': '',
                            'ko_number': annotation['ko_number'],
                            'kofamscan_annotation': annotation['definition']
                        })
        
        # Create summary DataFrame for this genome
        if summary_data:
            summary_df = pd.DataFrame(summary_data)
            
            # Sort by locus tag for better readability
            summary_df = summary_df.sort_values('locus_tag')
            
            # Save genome summary
            output_file = output_path / f"{genome}_genome_summary.tsv"
            summary_df.to_csv(output_file, sep='\t', index=False)
            print(f"  Saved {len(summary_df)} entries to {output_file}")
            print(f"  Added {annotations_added} annotations from result files")
            
            all_summaries.append(summary_df)
        else:
            print(f"  No data found for genome {genome}")
    
    # Create combined summary
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
        "--predicted-orfs-dir",
        required=True,
        help="Path to predicted_orfs directory containing .ffn files"
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
    
    # Validate predicted ORFs directory
    if not os.path.exists(args.predicted_orfs_dir):
        print(f"Error: Predicted ORFs directory {args.predicted_orfs_dir} does not exist!")
        sys.exit(1)
    
    # Create summaries
    create_genome_summaries(args.input_dir, args.output_dir, args.predicted_orfs_dir, args.genome_list)


if __name__ == "__main__":
    main() 