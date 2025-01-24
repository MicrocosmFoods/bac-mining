import argparse
from Bio import SeqIO
from collections import defaultdict
import os
import pandas as pd
import glob

def parse_args():
    parser = argparse.ArgumentParser(description="Split peptide predictions by genome and combine all predictions per genome.")
    
    parser.add_argument("--smorf-fasta", 
                        help="Combined FASTA file containing all smorf predictions")
    parser.add_argument("--encrypted-dir", 
                        help="Directory containing encrypted peptide prediction FASTA files")
    parser.add_argument("--cleavage-dir", 
                        help="Directory containing cleavage peptide prediction FASTA files")
    parser.add_argument("--ripp-dir", 
                        help="Directory containing RiPP core peptide prediction FASTA files")
    parser.add_argument("--genome-stb",
                        help="TSV file containing genome to scaffold mapping",
                        required=True)
    parser.add_argument("--outdir",
                        help="Output directory for per-genome FASTA files",
                        required=True)
    
    return parser.parse_args()

def read_fasta_file(fasta_file, tool_name):
    """Read single FASTA file and return dictionary of sequences with modified headers."""
    sequences = {}
    if fasta_file and os.path.exists(fasta_file):
        for record in SeqIO.parse(fasta_file, "fasta"):
            record.description=f"{record.description} tool={tool_name}"
            sequences[str(record.seq)] = record
    return sequences

def read_fasta_dir(directory, tool_name):
    """Read all FASTA files in a directory and return a dictionary of sequences with modified headers."""
    sequences = {}
    if directory and os.path.exists(directory):
        for fasta_file in glob.glob(os.path.join(directory, "*.fasta")):
            for record in SeqIO.parse(fasta_file, "fasta"):
                record.description=f"{record.description} tool={tool_name}"
                sequences[str(record.seq)] = record
    return sequences

def get_genome_name_from_header(header, genome_stb):
    """Extract genome name from header and use genome_stb to map to genome name."""



