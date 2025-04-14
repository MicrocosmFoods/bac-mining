#!/usr/bin/env python3

import os
import argparse
import math
from pathlib import Path

def create_genome_batches(genome_dir, batch_size, output_dir):
    # Get all genome files (supporting multiple extensions)
    genome_files = []
    for ext in ['.fa', '.fasta', '.fna']:
        genome_files.extend(list(Path(genome_dir).glob(f'*{ext}')))
    
    if not genome_files:
        raise ValueError(f"No genome files found in directory: {genome_dir}")
    
    # Sort the files to ensure consistent batching
    genome_files.sort()
    
    # Calculate number of batches needed
    total_genomes = len(genome_files)
    num_batches = math.ceil(total_genomes / batch_size)
    
    print(f"Found {total_genomes} genome files")
    print(f"Creating {num_batches} batches with approximately {batch_size} genomes each")
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Split genomes into batches and write to files
    for batch_num in range(num_batches):
        start_idx = batch_num * batch_size
        end_idx = min((batch_num + 1) * batch_size, total_genomes)
        
        batch_genomes = genome_files[start_idx:end_idx]
        output_file = os.path.join(output_dir, f"genomes_batch_{batch_num + 1}.txt")
        
        with open(output_file, 'w') as f:
            for genome in batch_genomes:
                # Write only the stem (basename without extension)
                f.write(f"{genome.stem}\n")
        
        print(f"Batch {batch_num + 1}: wrote {len(batch_genomes)} genomes to {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Split genome files into batches for processing")
    parser.add_argument("-i", "--input_dir", required=True,
                      help="Directory containing genome files (.fa, .fasta, or .fna)")
    parser.add_argument("-b", "--batch_size", type=int, required=True,
                      help="Number of genomes per batch")
    parser.add_argument("-o", "--output_dir", default="genome_batches",
                      help="Output directory for batch files (default: genome_batches)")
    
    args = parser.parse_args()
    
    create_genome_batches(args.input_dir, args.batch_size, args.output_dir)

if __name__ == "__main__":
    main()