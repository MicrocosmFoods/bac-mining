#!/usr/bin/env python

import argparse
import os
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description="Convert .gbk file to .gff format")
    parser.add_argument("--gbk", required=True, help="Input .gbk file")
    parser.add_argument("--output", required=True, help="Output .gff file")
    return parser.parse_args()

def gbk_to_gff(gbk_file, output_file):
    """Convert .gbk file to .gff format."""
    with open(output_file, 'w') as gff_out:
        # Write GFF header
        gff_out.write("##gff-version 3\n")
        
        for record in SeqIO.parse(gbk_file, "genbank"):
            seqid = record.id
            
            for feature in record.features:
                if feature.type == 'CDS':
                    # Get locus tag
                    locus_tag = feature.qualifiers.get('locus_tag', ['unknown'])[0]
                    
                    # Get coordinates
                    start = feature.location.start.position + 1  # GFF is 1-based
                    end = feature.location.end.position
                    
                    # Get strand
                    if feature.location.strand == 1:
                        strand = "+"
                    elif feature.location.strand == -1:
                        strand = "-"
                    else:
                        strand = "."
                    
                    # Get product annotation
                    product = feature.qualifiers.get('product', ['hypothetical protein'])[0]
                    
                    # Write GFF line
                    # Format: seqid source type start end score strand phase attributes
                    attributes = f"ID={locus_tag};locus_tag={locus_tag};product={product}"
                    gff_line = f"{seqid}\tpyrodigal\tCDS\t{start}\t{end}\t.\t{strand}\t0\t{attributes}\n"
                    gff_out.write(gff_line)

def main():
    args = parse_args()
    
    print(f"Converting {args.gbk} to {args.output}")
    gbk_to_gff(args.gbk, args.output)
    print(f"Conversion complete: {args.output}")

if __name__ == "__main__":
    main() 