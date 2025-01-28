#!/usr/bin/env python

import argparse
from Bio import SeqIO
from collections import defaultdict
import pandas as pd
import os
import glob

def parse_args():
    parser = argparse.ArgumentParser(description="Summarize peptide and BGC counts per genome")
    
    parser.add_argument("--peptide-dir",
                        help="Directory containing peptide FASTA files",
                        required=True)
    parser.add_argument("--bgc-summary",
                        help="TSV file containing BGC summary information",
                        required=True)
    parser.add_argument("--genome-stb",
                        help="TSV file mapping scaffolds to genome IDs",
                        required=True)
    parser.add_argument("--output",
                        help="Output TSV file path",
                        required=True)
    
    return parser.parse_args()

def count_peptide_types(fasta_file):
    """Count different peptide types in a FASTA file."""
    peptide_counts = defaultdict(int)
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        # Parse header which should be in format: genomename_peptidetype_XX
        parts = record.id.split('_')
        if len(parts) >= 2:
            peptide_type = parts[1]  # Get the peptide type
            peptide_counts[peptide_type] += 1
    
    return peptide_counts

def get_bgc_counts(bgc_summary_file, genome_stb_file):
    """Get BGC type counts from summary TSV using genome mapping."""
    # Read BGC summary and genome mapping files
    bgc_df = pd.read_csv(bgc_summary_file, sep='\t')
    genome_stb = pd.read_csv(genome_stb_file, sep='\t')
    
    # Merge BGC summary with genome mapping
    bgc_df = pd.merge(
        bgc_df,
        genome_stb[['scaffold_id', 'mag_id']],
        left_on='scaffold_name',  # column name in BGC summary
        right_on='scaffold_id',   # column name in genome STB
        how='left'
    )
    
    # Check for unmatched scaffolds
    unmatched = bgc_df[bgc_df['mag_id'].isna()]
    if not unmatched.empty:
        print("\nWarning: Some scaffolds in BGC summary could not be mapped to genomes:")
        for scaffold in unmatched['scaffold_name'].unique():
            print(f"  {scaffold}")
    
    # Split BGC types if they contain multiple types (e.g., "NRPS+T1PKS")
    bgc_df['bgc_type'] = bgc_df['bgc_type'].fillna('')
    
    # Create a dictionary to store counts per genome
    genome_bgc_counts = defaultdict(lambda: defaultdict(int))
    
    # Count each BGC type per genome
    for _, row in bgc_df.iterrows():
        if pd.notna(row['mag_id']) and row['bgc_type']:  # Skip entries without genome mapping or BGC type
            bgc_types = row['bgc_type'].split('+')
            for bgc_type in bgc_types:
                genome_bgc_counts[row['mag_id']][bgc_type.strip()] += 1
    
    return genome_bgc_counts

def main():
    args = parse_args()
    
    # Get all FASTA files in the directory
    fasta_files = glob.glob(os.path.join(args.peptide_dir, "*.fasta"))
    
    # Dictionary to store all counts per genome
    genome_counts = defaultdict(lambda: defaultdict(int))
    
    # Process each FASTA file
    for fasta_file in fasta_files:
        genome_name = os.path.basename(fasta_file).replace('.fasta', '')
        peptide_counts = count_peptide_types(fasta_file)
        
        # Store counts for this genome
        for peptide_type, count in peptide_counts.items():
            genome_counts[genome_name][peptide_type] = count
    
    # Get BGC counts using genome mapping
    bgc_counts = get_bgc_counts(args.bgc_summary, args.genome_stb)
    
    # Combine all counts into a DataFrame
    all_data = []
    
    # Get all possible column names
    peptide_types = set()
    bgc_types = set()
    
    # Collect all possible types
    for genome_data in genome_counts.values():
        peptide_types.update(genome_data.keys())
    for genome_data in bgc_counts.values():
        bgc_types.update(genome_data.keys())
    
    # Create sorted lists of types
    peptide_types = sorted(list(peptide_types))
    bgc_types = sorted(list(bgc_types))
    
    # Combine all data
    all_genomes = sorted(set(list(genome_counts.keys()) + list(bgc_counts.keys())))
    
    for genome in all_genomes:
        row_data = {'genome_name': genome}
        
        # Add peptide counts
        for ptype in peptide_types:
            row_data[ptype] = genome_counts[genome][ptype]
            
        # Add BGC counts
        for btype in bgc_types:
            row_data[f'bgc_{btype}'] = bgc_counts[genome][btype]
            
        all_data.append(row_data)
    
    # Create DataFrame and save to TSV
    df = pd.DataFrame(all_data)
    
    # Ensure genome_name is first column
    cols = ['genome_name'] + [col for col in df.columns if col != 'genome_name']
    df = df[cols]
    
    # Save to TSV
    df.to_csv(args.output, sep='\t', index=False)
    
    # Print summary
    print(f"\nSummary of counts:")
    print(f"Total genomes processed: {len(all_genomes)}")
    print("\nPeptide types found:")
    for ptype in peptide_types:
        total = df[ptype].sum()
        genomes_with_type = (df[ptype] > 0).sum()
        print(f"  {ptype}: {total} total peptides in {genomes_with_type} genomes")
    
    print("\nBGC types found:")
    for btype in bgc_types:
        col = f'bgc_{btype}'
        total = df[col].sum()
        genomes_with_type = (df[col] > 0).sum()
        print(f"  {btype}: {total} total BGCs in {genomes_with_type} genomes")

if __name__ == "__main__":
    main()