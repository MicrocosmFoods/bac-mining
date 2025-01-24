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

def get_genome_from_header(header, genome_mapping):
    """Extract genome ID from sequence header using genome mapping."""
    header_parts = header.split()[0]
    
    for genome_id in genome_mapping:
        if genome_id in header_parts:
            return genome_id
    return None

def main():
    args = parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.outdir, exist_ok=True)
    
    # Read genome mapping
    genome_stb = pd.read_csv(args.genome_stb, sep='\t')
    genome_mapping = set(genome_stb['mag_id'].unique())
    
    # Dictionary to store sequences by genome
    genome_sequences = defaultdict(dict)
    
    # Process smorf proteins (single file)
    smorf_sequences = read_fasta_file(args.smorf_fasta, 'smorf')
    for seq, record in smorf_sequences.items():
        genome_id = get_genome_from_header(record.description, genome_mapping)
        if genome_id and seq not in genome_sequences[genome_id]:
            genome_sequences[genome_id][seq] = record
    
    # Process other peptide types (directories)
    dir_mapping = {
        'encrypted': args.encrypted_dir,
        'cleavage': args.cleavage_dir,
        'ripp': args.ripp_dir
    }
    
    for tool_name, directory in dir_mapping.items():
        sequences = read_fasta_dir(directory, tool_name)
        for seq, record in sequences.items():
            genome_id = get_genome_from_header(record.description, genome_mapping)
            if genome_id and seq not in genome_sequences[genome_id]:
                genome_sequences[genome_id][seq] = record

    # Write combined FASTA files for each genome
    for genome_id, sequences in genome_sequences.items():
        output_file = os.path.join(args.outdir, f"{genome_id}_combined_peptides.fasta")
        with open(output_file, 'w') as f:
            SeqIO.write(sequences.values(), f, "fasta")
        print(f"Written {len(sequences)} unique peptides for genome {genome_id}")

if __name__ == "__main__":
    main()
