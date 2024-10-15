#!/usr/bin/env python3

import argparse
from collections import defaultdict
from Bio import SeqIO

def count_proteins(fasta_file):
    """
    Count the number of proteins per genome in a FASTA file.
    """
    counts = defaultdict(int)
    for record in SeqIO.parse(fasta_file, "fasta"):
        genome_name = record.id.split("_id_")[0]
        counts[genome_name] += 1
    return counts

def parse_args():
    parser = argparse.ArgumentParser(description="Count smORF proteins per genome in combined and non-redundant FASTA files.")
    parser.add_argument("combined_fasta", help="Path to the combined FASTA file")
    parser.add_argument("nonredundant_fasta", help="Path to the non-redundant FASTA file")
    parser.add_argument("output_tsv", help="Path to the output TSV file")
    return parser.parse_args()

def main():
    args = parse_args()

    # Count proteins in combined FASTA
    combined_counts = count_proteins(args.combined_fasta)

    # Count proteins in non-redundant FASTA
    nonredundant_counts = count_proteins(args.nonredundant_fasta)

    # Combine results and write to TSV
    with open(args.output_tsv, 'w') as outfile:
        outfile.write("Genome\tCombined_smORFs\tNonredundant_smORFs\n")
        for genome in sorted(set(combined_counts.keys()) | set(nonredundant_counts.keys())):
            combined = combined_counts.get(genome, 0)
            nonredundant = nonredundant_counts.get(genome, 0)
            outfile.write(f"{genome}\t{combined}\t{nonredundant}\n")

    print(f"Results written to {args.output_tsv}")

if __name__ == "__main__":
    main()