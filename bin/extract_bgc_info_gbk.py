#!/usr/bin/env python

import argparse
import os
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

class AntismashGBKParser:
    """Parse GBK file from antiSMASH to extract relevant BGC and peptide information."""

    def __init__(self, gbk_files):
        self.gbk_files = gbk_files

    @property
    def bgc_length(self):
        """Return the total length of all BGCs."""
        total_length = 0
        for gbk_file in self.gbk_files:
            for seq_record in SeqIO.parse(gbk_file, 'genbank'):
                total_length += len(seq_record.seq)
        return total_length

    @property
    def bgc_type(self):
        """Return the types of BGC regions."""
        bgc_type_list = []
        for gbk_file in self.gbk_files:
            for seq_record in SeqIO.parse(gbk_file, 'genbank'):
                for feature in seq_record.features:
                    if feature.type == 'region':
                        bgc_type_list += feature.qualifiers.get('product', [])
        return '+'.join(set(bgc_type_list))

    @property
    def bgc_completeness(self):
        """Return the completeness of the BGCs."""
        for gbk_file in self.gbk_files:
            for seq_record in SeqIO.parse(gbk_file, 'genbank'):
                for feature in seq_record.features:
                    if feature.type == 'region':
                        contig_edge = feature.qualifiers.get("contig_edge", ["complete"])[0]
                        if contig_edge != "complete":
                            return "incomplete"
        return "complete"

    def ex_lanthipeptide(self):
        """Extract leader and core sequences of lanthipeptides (RiPPs)."""
        leader_seq, core_seq = [], []
        for gbk_file in self.gbk_files:
            for seq_record in SeqIO.parse(gbk_file, 'genbank'):
                for feature in seq_record.features:
                    if feature.type == 'CDS_motif':
                        qualifiers = feature.qualifiers
                        if qualifiers.get("core_sequence") and qualifiers.get("leader_sequence"):
                            leader_seq += qualifiers.get("leader_sequence")
                            core_seq += qualifiers.get("core_sequence")
        return leader_seq, core_seq

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Parse antiSMASH GBK files to summarize BGCs and extract lanthipeptides.")
    parser.add_argument("genome_name", help="Name of the genome (input argument)")
    parser.add_argument("input_gbk", nargs='+', help="Paths to the input GBK files")
    parser.add_argument("output_bgc_tsv", help="Path to the output TSV file for BGC information")
    parser.add_argument("output_peptide_tsv", help="Path to the output TSV file for lanthipeptide information")
    parser.add_argument("output_peptide_fasta", help="Path to the output FASTA file for core peptide sequences")
    return parser.parse_args()

def main():
    args = parse_args()

    # Parse the GBK files
    parser = AntismashGBKParser(args.input_gbk)
    bgc_type = parser.bgc_type
    bgc_length = parser.bgc_length
    bgc_completeness = parser.bgc_completeness

    # Extract lanthipeptides (RiPPs)
    leader_seqs, core_seqs = parser.ex_lanthipeptide()

    # Write the BGC data to the first TSV file
    bgc_data = {
        'genome_name': args.genome_name,
        'bgc_type': bgc_type,
        'bgc_length': bgc_length,
        'bgc_completeness': bgc_completeness,
    }
    bgc_df = pd.DataFrame([bgc_data])
    bgc_df.to_csv(args.output_bgc_tsv, sep='\t', index=False)

    # If lanthipeptides are found, write them to the second TSV file and FASTA file
    if leader_seqs and core_seqs:
        peptide_data = {
            'genome_name': args.genome_name,
            'leader_seq': ';'.join(leader_seqs),
            'core_seq': ';'.join(core_seqs),
        }
        peptide_df = pd.DataFrame([peptide_data])
        peptide_df.to_csv(args.output_peptide_tsv, sep='\t', index=False)
        print(f"Lanthipeptides found and written to {args.output_peptide_tsv}.")

        # Create and write peptide FASTA file
        fasta_records = []
        for i, seq in enumerate(core_seqs, 1):
            record = SeqRecord(Seq(seq), 
                               id=f"{args.genome_name}_{i}", 
                               description="")
            fasta_records.append(record)
        
        SeqIO.write(fasta_records, args.output_peptide_fasta, "fasta")
        print(f"Core peptide sequences written to {args.output_peptide_fasta}.")
    else:
        print(f"No lanthipeptides found in the input GBK files. No peptide files generated.")

if __name__ == "__main__":
    main()