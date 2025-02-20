#!/usr/bin/env python3

import argparse
from Bio import SeqIO

def filter_small_proteins(input_fasta, output_fasta, max_length=50):
    """Filter proteins to keep only those shorter than max_length."""
    with open(input_fasta, 'r') as in_handle, open(output_fasta, 'w') as out_handle:
        for record in SeqIO.parse(in_handle, 'fasta'):
            if len(record.seq) <= max_length:
                SeqIO.write(record, out_handle, 'fasta')

def main():
    parser = argparse.ArgumentParser(description='Filter proteins by length')
    parser.add_argument('input_fasta', help='Input FASTA file')
    parser.add_argument('output_fasta', help='Output FASTA file')
    parser.add_argument('--max_length', type=int, default=50,
                       help='Maximum protein length to keep (default: 50)')
    args = parser.parse_args()

    filter_small_proteins(args.input_fasta, args.output_fasta, args.max_length)

if __name__ == '__main__':
    main()