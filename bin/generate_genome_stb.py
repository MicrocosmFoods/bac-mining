import os
import argparse
from Bio import SeqIO
from collections import defaultdict

def parse_fasta_files(fasta_files, output_file):
    # Track scaffold IDs and their corresponding MAG IDs
    scaffold_to_mags = defaultdict(list)
    # Track scaffold IDs within each MAG to check for internal duplicates
    mag_scaffolds = defaultdict(set)
    
    # First pass: collect all scaffold IDs and their MAGs
    for fasta_file in fasta_files:
        mag_id = os.path.splitext(os.path.basename(fasta_file))[0]
        with open(fasta_file, 'r') as f:
            for record in SeqIO.parse(f, "fasta"):
                scaffold_id = record.id
                scaffold_to_mags[scaffold_id].append(mag_id)
                
                # Check for duplicates within the same MAG
                if scaffold_id in mag_scaffolds[mag_id]:
                    duplicates = {scaffold_id: [mag_id]}
                    duplicate_file = output_file + ".duplicates.txt"
                    with open(duplicate_file, 'w') as dup_f:
                        dup_f.write("scaffold_id\tmag_ids\tnote\n")
                        dup_f.write(f"{scaffold_id}\t{mag_id}\tduplicate within same MAG\n")
                    raise AssertionError(f"Duplicate scaffold ID '{scaffold_id}' found within MAG '{mag_id}'. Details written to: {duplicate_file}")
                
                mag_scaffolds[mag_id].add(scaffold_id)
    
    # Check for duplicate scaffold IDs across different MAGs
    duplicates = {scaffold: mags for scaffold, mags in scaffold_to_mags.items() if len(mags) > 1}
    
    if duplicates:
        # Write duplicates to a file
        duplicate_file = output_file + ".duplicates.txt"
        with open(duplicate_file, 'w') as dup_f:
            dup_f.write("scaffold_id\tmag_ids\tnote\n")
            for scaffold, mags in duplicates.items():
                dup_f.write(f"{scaffold}\t{', '.join(mags)}\tduplicate across MAGs\n")
        
        error_msg = f"Duplicate scaffold IDs found across different MAGs. Details written to: {duplicate_file}"
        raise AssertionError(error_msg)
    
    # If no duplicates, write the output file
    with open(output_file, 'w') as out_tsv:
        out_tsv.write("mag_id\tscaffold_id\n")
        for fasta_file in fasta_files:
            mag_id = os.path.splitext(os.path.basename(fasta_file))[0]
            with open(fasta_file, 'r') as f:
                for record in SeqIO.parse(f, "fasta"):
                    scaffold_id = record.id
                    out_tsv.write(f"{mag_id}\t{scaffold_id}\n")

def main():
    parser = argparse.ArgumentParser(description="Parse multiple FASTA files and create a TSV with filenames and scaffold names.")
    parser.add_argument("fasta_files", nargs='+', help="List of FASTA files to process")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file")

    args = parser.parse_args()

    parse_fasta_files(args.fasta_files, args.output)

if __name__ == "__main__":
    main()