#!/usr/bin/env python3

from Bio import SeqIO
import os
import argparse

def get_genome_name(filename):
    return os.path.splitext(os.path.basename(filename))[0]


def process_fasta_file(input_file, genome_name):
    new_records = {}
    with open(input_file, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            # Remove trailing '*' if present
            if record.seq.endswith('*'):
                record.seq = record.seq[:-1]
            record.id = f"{genome_name}_id_{record.id}"
            new_records[record.id] = record
    return new_records


def combine_fastas(input_files, output_file):
    all_records = {}
    for input_file in input_files:
        genome_name = get_genome_name(input_file)
        new_records = process_fasta_file(input_file, genome_name)
        all_records.update(new_records)

    with open(output_file, "w") as outfile:
        SeqIO.write(all_records.values(), outfile, "fasta")

def parse_arguments():
    parser = argparse.ArgumentParser(description='Rename FASTA headers for input files and concatenate to single FASTA')
    parser.add_argument('input_files', nargs='+', help='Input FASTA files')
    parser.add_argument('output_file', help='Output concatenated FASTA file')
    return parser.parse_args()

def main():
    args = parse_arguments()
    combine_fastas(args.input_files, args.output_file)


if __name__ == '__main__':
    main()