#!/usr/bin/env python

import argparse
from Bio import SeqIO
from collections import defaultdict
import os
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description="Split peptide predictions by genome and combine all predictions per genome.")
    
    parser.add_argument("--smorf-fasta", 
                        help="Combined FASTA file containing all smorf predictions")
    parser.add_argument("--encrypted-fasta", 
                        help="Combined FASTA file containing all encrypted peptide predictions")
    parser.add_argument("--cleavage-fasta", 
                        help="Combined FASTA file containing all cleavage peptide predictions")
    parser.add_argument("--ripp-fasta", 
                        help="Combined FASTA file containing all RiPP core peptide predictions")
    parser.add_argument("--genome-stb",
                        help="TSV file containing genome to scaffold mapping",
                        required=True)
    parser.add_argument("--outdir",
                        help="Output directory for per-genome FASTA files",
                        required=True)
    
    return parser.parse_args()

def read_fasta_file(fasta_file, tool_name, scaffold_to_genome=None):
    """Read FASTA file and return dictionary of sequences grouped by genome."""
    sequences = defaultdict(list)
    if fasta_file and os.path.exists(fasta_file):
        for record in SeqIO.parse(fasta_file, "fasta"):
            if tool_name == 'ripp' and scaffold_to_genome is not None:
                # For RiPPs, get scaffold ID from header and map to genome
                scaffold_id = record.id.split('_')[0]
                if scaffold_id in scaffold_to_genome:
                    genome_id = scaffold_to_genome[scaffold_id]
                    sequences[genome_id].append((str(record.seq), record))
            else:
                # For other peptides, genome ID is already in header
                genome_id = record.description.split()[0]
                sequences[genome_id].append((str(record.seq), record))
    return sequences

def create_new_record(seq, genome_id, tool_name, count):
    """Create a new SeqRecord with formatted header."""
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    
    new_id = f"{genome_id}_{tool_name}_peptide_{count:03d}"
    return SeqRecord(
        Seq(seq),
        id=new_id,
        name=new_id,
        description=""
    )

def main():
    args = parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.outdir, exist_ok=True)
    
    # Read genome mapping
    genome_stb = pd.read_csv(args.genome_stb, sep='\t')
    scaffold_to_genome = dict(zip(genome_stb['scaffold_id'], genome_stb['mag_id']))
    
    # Dictionary to store sequences by genome
    genome_sequences = defaultdict(lambda: defaultdict(list))
    
    # Process each peptide type
    peptide_files = {
        'smorf': args.smorf_fasta,
        'encrypted': args.encrypted_fasta,
        'cleavage': args.cleavage_fasta,
        'ripp': args.ripp_fasta
    }
    
    # Read all sequences and group by genome and tool
    for tool_name, fasta_file in peptide_files.items():
        scaffold_map = scaffold_to_genome if tool_name == 'ripp' else None
        sequences = read_fasta_file(fasta_file, tool_name, scaffold_map)
        
        # Group sequences by genome
        for genome_id, seqs in sequences.items():
            # Add only unique sequences for this tool and genome
            seen_seqs = set()
            for seq, record in seqs:
                if seq not in seen_seqs:
                    genome_sequences[genome_id][tool_name].append(seq)
                    seen_seqs.add(seq)
    
    # Write combined FASTA files for each genome
    for genome_id, tool_sequences in genome_sequences.items():
        output_file = os.path.join(args.outdir, f"{genome_id}.fasta")
        records = []
        
        # Process each tool's sequences
        for tool_name, seqs in tool_sequences.items():
            # Create new records with formatted headers
            for i, seq in enumerate(seqs, 1):
                record = create_new_record(seq, genome_id, tool_name, i)
                records.append(record)
        
        # Write all records for this genome
        with open(output_file, 'w') as f:
            SeqIO.write(records, f, "fasta")
        
        # Print summary
        tool_counts = {tool: len(seqs) for tool, seqs in tool_sequences.items()}
        print(f"Genome {genome_id}:")
        for tool, count in tool_counts.items():
            print(f"  {tool}: {count} peptides")
        print(f"  Total: {sum(tool_counts.values())} peptides")
        print("---")

if __name__ == "__main__":
    main()