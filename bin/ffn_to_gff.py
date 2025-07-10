#!/usr/bin/env python

import argparse
import os
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description="Convert .ffn file to .gff format")
    parser.add_argument("--ffn", required=True, help="Input .ffn file")
    parser.add_argument("--output", required=True, help="Output .gff file")
    return parser.parse_args()

def extract_location_from_header(description):
    """Extract location information from FASTA header description."""
    if "[location=" in description:
        location_part = description.split("[location=")[1].split("]")[0]
        return location_part
    return None

def parse_location(location_str):
    """Parse location string to get start, end, and strand."""
    if not location_str:
        return None, None, None
    
    # Handle complement locations
    strand = "+"
    if "complement(" in location_str:
        strand = "-"
        location_str = location_str.replace("complement(", "").replace(")", "")
    
    # Handle join locations (multiple segments)
    if "join(" in location_str:
        # For join locations, use the first segment
        location_str = location_str.replace("join(", "").replace(")", "")
        segments = location_str.split(",")
        if segments:
            location_str = segments[0]
    
    # Parse start..end
    if ".." in location_str:
        start_str, end_str = location_str.split("..")
        try:
            start = int(start_str)
            end = int(end_str)
            return start, end, strand
        except ValueError:
            return None, None, None
    
    return None, None, None

def ffn_to_gff(ffn_file, output_file):
    """Convert .ffn file to .gff format."""
    with open(output_file, 'w') as gff_out:
        # Write GFF header
        gff_out.write("##gff-version 3\n")
        
        for record in SeqIO.parse(ffn_file, "fasta"):
            locus_tag = record.id
            description = record.description
            
            # Extract location from header
            location_str = extract_location_from_header(description)
            start, end, strand = parse_location(location_str)
            
            if start is not None and end is not None:
                # Write GFF line
                # Format: seqid source type start end score strand phase attributes
                gff_line = f"{locus_tag}\tpyrodigal\tCDS\t{start}\t{end}\t.\t{strand}\t0\tID={locus_tag};locus_tag={locus_tag}\n"
                gff_out.write(gff_line)

def main():
    args = parse_args()
    
    print(f"Converting {args.ffn} to {args.output}")
    ffn_to_gff(args.ffn, args.output)
    print(f"Conversion complete: {args.output}")

if __name__ == "__main__":
    main() 