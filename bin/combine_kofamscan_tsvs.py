#!/usr/bin/env python

import argparse
import polars as pl
import os

def parse_kofamscan_file(file_path):
    """Parse a single KofamScan TSV file and return significant hits with genome name."""
    # Get genome name from file name
    genome_name = os.path.basename(file_path).replace('_kofamscan_annotations.tsv', '')
    
    # Read file, keeping only significant hits
    significant_hits = []
    with open(file_path, 'r') as f:
        # Skip header lines
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('*'):
                # Remove the * and split the line
                fields = line.strip()[1:].split(maxsplit=5)  # Split only first 5 fields
                if len(fields) == 6:  # Ensure we have all fields
                    significant_hits.append(fields)
    
    # Create dataframe
    if significant_hits:
        df = pl.DataFrame(
            significant_hits,
            schema=['gene_name', 'KO', 'threshold', 'score', 'E-value', 'KO_definition']
        )
        # Add genome name
        df = df.with_columns(pl.lit(genome_name).alias('genome_name'))
        return df
    return None

def main():
    parser = argparse.ArgumentParser(description="Combine KofamScan results from multiple genomes")
    parser.add_argument(
        "--input_files",
        nargs='+',
        required=True,
        help="Input KofamScan TSV files"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output combined TSV file"
    )
    
    args = parser.parse_args()
    
    # Process all input files
    all_dfs = []
    for file_path in args.input_files:
        df = parse_kofamscan_file(file_path)
        if df is not None:
            all_dfs.append(df)
    
    # Combine all dataframes
    if all_dfs:
        combined_df = pl.concat(all_dfs)
        # Reorder columns to put genome_name first
        cols = ['genome_name'] + [col for col in combined_df.columns if col != 'genome_name']
        combined_df = combined_df.select(cols)
        # Write output
        combined_df.write_csv(args.output, separator='\t')
        
        # Print summary
        print(f"\nProcessed {len(args.input_files)} input files")
        
if __name__ == "__main__":
    main()