#!/usr/bin/env python

import argparse
import pandas as pd
import os
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description="Summarize peptide and BGC counts per genome")
    
    parser.add_argument("--smorfinder-tsvs",
                        nargs='+',
                        help="List of smorfinder TSV files",
                        required=True)
    parser.add_argument("--deeppeptide-tsvs",
                        nargs='+',
                        help="List of deeppeptide TSV files",
                        required=True)
    parser.add_argument("--bgc-summary",
                        help="TSV file containing BGC summary information",
                        required=True)
    parser.add_argument("--output-counts",
                        help="Output TSV file path for molecule counts",
                        required=True)
    parser.add_argument("--output-smorfinder",
                        help="Output TSV file path for combined smorfinder results",
                        required=True)
    parser.add_argument("--output-deeppeptide",
                        help="Output TSV file path for combined deeppeptide results",
                        required=True)
    
    return parser.parse_args()

def combine_prediction_tsvs(tsv_files, output_file):
    """Combine prediction TSVs and add genome names."""
    all_dfs = []
    for tsv_file in tsv_files:
        # Remove both _smorfinder.tsv and _deeppeptide.tsv endings
        genome_name = os.path.basename(tsv_file)
        genome_name = genome_name.replace('_smorfinder.tsv', '').replace('_deeppeptide.tsv', '')
        df = pd.read_csv(tsv_file, sep='\t')
        df['genome_name'] = genome_name
        all_dfs.append(df)
    
    if all_dfs:
        combined_df = pd.concat(all_dfs, ignore_index=True)
        combined_df.to_csv(output_file, sep='\t', index=False)

def classify_bgc_type(bgc_type):
    """Classify BGC type into main categories or 'other'."""
    main_types = {'NRPS', 'T1PKS', 'T3PKS', 'RiPP-like', 'betalactone', 'terpene'}
    
    # If it contains a '+', it's a hybrid - classify as other
    if '+' in bgc_type:
        return 'other'
    
    # If it's one of our main types, return it
    if bgc_type.strip() in main_types:
        return bgc_type.strip()
    
    # Everything else is other
    return 'other'

def count_molecules_per_genome():
    """Count different molecule types per genome from TSV files."""
    args = parse_args()
    
    # Read BGC summary
    bgc_df = pd.read_csv(args.bgc_summary, sep='\t')
    
    # Count BGCs per genome with simplified categories
    # Group by mag_id and bgc_file to count unique BGCs, not genes
    bgc_counts = defaultdict(lambda: defaultdict(int))
    unique_bgcs = set()
    
    for _, row in bgc_df.iterrows():
        if pd.notna(row['bgc_type']) and pd.notna(row['mag_id']) and pd.notna(row['bgc_file']):
            # Create unique BGC identifier
            bgc_id = f"{row['mag_id']}_{row['bgc_file']}"
            
            # Only count each unique BGC once
            if bgc_id not in unique_bgcs:
                unique_bgcs.add(bgc_id)
                bgc_type = classify_bgc_type(row['bgc_type'])
                bgc_counts[row['mag_id']][bgc_type] += 1
    
    # Process smorfinder results
    smorf_counts = defaultdict(lambda: defaultdict(int))
    for tsv_file in args.smorfinder_tsvs:
        genome_name = os.path.basename(tsv_file).replace('_smorfinder.tsv', '')
        df = pd.read_csv(tsv_file, sep='\t')
        smorf_counts[genome_name]['smorf'] = len(df)
    
    # Process deeppeptide results
    deep_counts = defaultdict(lambda: defaultdict(int))
    for tsv_file in args.deeppeptide_tsvs:
        genome_name = os.path.basename(tsv_file).replace('_deeppeptide.tsv', '')
        df = pd.read_csv(tsv_file, sep='\t')
        for _, row in df.iterrows():
            deep_counts[genome_name][row['peptide_class']] += 1
    
    # Combine all counts
    all_genomes = set()
    all_genomes.update(bgc_counts.keys(), smorf_counts.keys(), deep_counts.keys())
    
    # Get all column types
    bgc_types = {'NRPS', 'T1PKS', 'T3PKS', 'RiPP-like', 'betalactone', 'bacteriocin', 'terpene', 'other'}
    deep_types = set()
    for genome_counts in deep_counts.values():
        deep_types.update(genome_counts.keys())
    
    # Create final dataframe
    rows = []
    for genome in sorted(all_genomes):
        row = {'genome_name': genome}
        
        # Add BGC counts
        for btype in sorted(bgc_types):
            row[f'bgc_{btype}'] = bgc_counts[genome][btype]
        
        # Add smorfinder count
        row['smorf'] = smorf_counts[genome]['smorf']
        
        # Add deeppeptide counts
        for dtype in sorted(deep_types):
            row[f'deeppeptide_{dtype}'] = deep_counts[genome][dtype]
        
        rows.append(row)
    
    # Create and save counts dataframe
    counts_df = pd.DataFrame(rows)
    counts_df.to_csv(args.output_counts, sep='\t', index=False)
    
    # Combine prediction TSVs
    combine_prediction_tsvs(args.smorfinder_tsvs, args.output_smorfinder)
    combine_prediction_tsvs(args.deeppeptide_tsvs, args.output_deeppeptide)

if __name__ == "__main__":
    count_molecules_per_genome()