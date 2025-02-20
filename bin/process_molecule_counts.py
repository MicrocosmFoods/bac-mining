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
    parser.add_argument("--antismash-peptides",
                        help="TSV file containing antismash peptide core sequences",
                        required=False)
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
    main_types = {'NRPS', 'T1PKS', 'T3PKS', 'RiPP-like', 'betalactone', 'bacteriocin', 'terpene'}
    
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

    # Only process antismash peptides if file is provided
    antismash_counts = defaultdict(int)
    if args.antismash_peptides and os.path.exists(args.antismash_peptides):
        peptides_df = pd.read_csv(args.antismash_peptides, sep='\t')
        for _, row in peptides_df.iterrows():
            genome_name = row['genome_name']
            antismash_counts[genome_name] += 1

    # Create summary dataframe
    summary_data = []
    for genome in set(list(smorf_counts.keys()) + list(deep_counts.keys()) + list(antismash_counts.keys())):
        row = {
            'genome': genome,
            'total_smorfs': smorf_counts[genome]['smorf'],
            'total_deeppeptide': sum(deep_counts[genome].values()),
        }
        
        # Add BGC counts
        genome_bgcs = bgc_df[bgc_df['genome_name'] == genome]
        row['total_bgcs'] = len(genome_bgcs) if not genome_bgcs.empty else 0
        
        # Only add antismash peptides if we processed them
        if args.antismash_peptides:
            row['total_antismash_peptides'] = antismash_counts[genome]
        
        summary_data.append(row)

    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(args.output_counts, sep='\t', index=False)
    
    # Combine prediction TSVs
    combine_prediction_tsvs(args.smorfinder_tsvs, args.output_smorfinder)
    combine_prediction_tsvs(args.deeppeptide_tsvs, args.output_deeppeptide)

if __name__ == "__main__":
    count_molecules_per_genome()