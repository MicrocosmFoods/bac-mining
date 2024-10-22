import os
import argparse
from Bio import SeqIO

def parse_fasta_files(fasta_files, output_file):
    with open(output_file, 'w') as out_tsv:
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