#!/usr/bin/env python3

import os
import argparse
from Bio import SeqIO
from collections import defaultdict

def check_and_fix_duplicates(fasta_dir, output_dir, duplicates_file):
    # Get all FASTA files from the directory
    fasta_files = [os.path.join(fasta_dir, f) for f in os.listdir(fasta_dir) 
                   if f.endswith(('.fa', '.fasta', '.fna'))]
    
    if not fasta_files:
        raise ValueError(f"No FASTA files found in directory: {fasta_dir}")
    
    # Track scaffold IDs and their corresponding MAG IDs
    scaffold_to_mags = defaultdict(list)
    # Track scaffold IDs within each MAG to check for internal duplicates
    mag_scaffolds = defaultdict(set)
    # Store records that need to be modified
    records_to_fix = []
    
    # collect all scaffold IDs and their MAGs
    for fasta_file in fasta_files:
        mag_id = os.path.splitext(os.path.basename(fasta_file))[0]
        with open(fasta_file, 'r') as f:
            for record in SeqIO.parse(f, "fasta"):
                scaffold_id = record.id
                scaffold_to_mags[scaffold_id].append((mag_id, record, fasta_file))
                
                # Check for duplicates within the same MAG
                if scaffold_id in mag_scaffolds[mag_id]:
                    records_to_fix.append((scaffold_id, mag_id, record, fasta_file, "within_mag"))
                
                mag_scaffolds[mag_id].add(scaffold_id)
    
    # Check for duplicate scaffold IDs across different MAGs
    for scaffold_id, mag_records in scaffold_to_mags.items():
        if len(mag_records) > 1:
            for mag_id, record, fasta_file in mag_records:
                records_to_fix.append((scaffold_id, mag_id, record, fasta_file, "across_mags"))
    
    # if found duplicates, write them to file and create fixed FASTA files
    if records_to_fix:
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Write duplicates report
        with open(duplicates_file, 'w') as dup_f:
            dup_f.write("mag_id\toriginal_scaffold_id\tnew_scaffold_id\n")
            
            # Track which files we need to rewrite
            files_to_rewrite = defaultdict(dict)
            
            for scaffold_id, mag_id, record, fasta_file, dup_type in records_to_fix:
                # Create new scaffold ID by appending MAG name
                new_scaffold_id = f"{scaffold_id}_{mag_id}"
                
                # Write to duplicates report
                dup_f.write(f"{mag_id}\t{scaffold_id}\t{new_scaffold_id}\n")
                
                # Store the modified record for rewriting, keyed by original scaffold ID
                record.id = new_scaffold_id
                record.description = ""  # Clear description to avoid duplicate IDs in header
                files_to_rewrite[fasta_file][scaffold_id] = record
            
            # Rewrite FASTA files that had duplicates
            for fasta_file in files_to_rewrite:
                # Read all records from original file
                output_records = []
                with open(fasta_file, 'r') as f:
                    for record in SeqIO.parse(f, "fasta"):
                        # If this is a scaffold that needs to be renamed, use the new version
                        if record.id in files_to_rewrite[fasta_file]:
                            output_records.append(files_to_rewrite[fasta_file][record.id])
                        else:
                            output_records.append(record)
                
                # Write new file with fixed records
                output_file = os.path.join(output_dir, os.path.basename(fasta_file))
                with open(output_file, 'w') as f:
                    SeqIO.write(output_records, f, "fasta")
        
        print(f"Found duplicate scaffold IDs. Details written to: {duplicates_file}")
        print(f"Fixed FASTA files written to: {output_dir}")
    else:
        print("No duplicates found in the FASTA files.")

def main():
    parser = argparse.ArgumentParser(description="Check for duplicate scaffold IDs in FASTA files and create fixed versions.")
    parser.add_argument("-i", "--input_dir", required=True, help="Directory containing FASTA files")
    parser.add_argument("-o", "--output_dir", required=True, help="Output directory for fixed FASTA files")
    parser.add_argument("-d", "--duplicates", required=True, help="Output file for duplicate entries")

    args = parser.parse_args()
    check_and_fix_duplicates(args.input_dir, args.output_dir, args.duplicates)

if __name__ == "__main__":
    main()