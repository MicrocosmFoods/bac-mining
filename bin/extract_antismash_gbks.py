#!/usr/bin/env python

import argparse
import os
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

class AntismashGBKParser:
    """Parse GBK file from antiSMASH to extract relevant BGC and peptide information."""

    def __init__(self, gbk_file, mag_scaffold_df):
        self.gbk_file = gbk_file
        self.records = list(SeqIO.parse(self.gbk_file, 'genbank'))
        self.mag_scaffold_df = mag_scaffold_df
        
    def parse_records(self):
        """Parse all records in the GBK file."""
        bgc_data_list = []
        peptide_data_list = []
        fasta_records = []
        for record in self.records:
            scaffold_name = record.id
            bgc_type = self.get_bgc_type(record)
            bgc_length = len(record.seq)
            leader_seqs, core_seqs = self.ex_lanthipeptide(record)
            bgc_data_list.append({
                'scaffold_name': scaffold_name,
                'bgc_file': os.path.basename(self.gbk_file),
                'bgc_type': bgc_type,
                'bgc_length': bgc_length,
            })
            if leader_seqs and core_seqs:
                peptide_data_list.append({
                    'scaffold_name': scaffold_name,
                    'bgc_file': os.path.basename(self.gbk_file),
                    'leader_seq': ';'.join(leader_seqs),
                    'core_seq': ';'.join(core_seqs),
                })
                
                # Look up the genome name from the mapping dataframe
                genome_name = self.mag_scaffold_df[
                    self.mag_scaffold_df['scaffold_id'] == scaffold_name
                ]['mag_id'].iloc[0] if not self.mag_scaffold_df[
                    self.mag_scaffold_df['scaffold_id'] == scaffold_name
                ].empty else 'unknown_genome'
                
                for i, seq in enumerate(core_seqs, 1):
                    record = SeqRecord(
                        Seq(seq),
                        id=f"{genome_name}_{scaffold_name}_{i}",
                        description=""
                    )
                    fasta_records.append(record)
        return bgc_data_list, peptide_data_list, fasta_records
    
    def get_bgc_type(self, record):
        """Return the type of BGC region."""
        bgc_type_list = []
        for feature in record.features:
            if feature.type == 'region':
                bgc_type_list += feature.qualifiers.get('product', [])
        return '+'.join(set(bgc_type_list))

    def ex_lanthipeptide(self, record):
        """Extract leader and core sequences of lanthipeptides (RiPPs)."""
        leader_seq, core_seq = [], []
        for feature in record.features:
            if feature.type == 'CDS_motif':
                qualifiers = feature.qualifiers
                if qualifiers.get("core_sequence") and qualifiers.get("leader_sequence"):
                    leader_seq += qualifiers.get("leader_sequence")
                    core_seq += qualifiers.get("core_sequence")
        return leader_seq, core_seq

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Parse antiSMASH GBK files to summarize BGCs and extract lanthipeptides.")
    parser.add_argument("input_gbk", nargs='+', help="Paths to the input GBK files")
    parser.add_argument("mag_scaffold_tsv", help="Path to the TSV file containing mag_id and scaffold_id mapping")
    parser.add_argument("output_bgc_tsv", help="Path to the output TSV file for BGC information")
    parser.add_argument("output_peptide_tsv", help="Path to the output TSV file for lanthipeptide information")
    parser.add_argument("output_peptide_fasta", help="Path to the output FASTA file for lanthipeptide sequences")
    return parser.parse_args()

def main():
    args = parse_args()

    # Read the mag_scaffold mapping file
    mag_scaffold_df = pd.read_csv(args.mag_scaffold_tsv, sep='\t')

    all_bgc_data = []
    all_peptide_data = []
    all_fasta_records = []
    for gbk_file in args.input_gbk:
        parser = AntismashGBKParser(gbk_file, mag_scaffold_df)
        bgc_data, peptide_data, fasta_records = parser.parse_records()
        all_bgc_data.extend(bgc_data)
        all_peptide_data.extend(peptide_data)
        all_fasta_records.extend(fasta_records)

    # Write BGC data to TSV
    bgc_df = pd.DataFrame(all_bgc_data)
    
    # Merge with mag_scaffold mapping
    bgc_df = pd.merge(bgc_df, mag_scaffold_df, 
                     left_on='scaffold_name', 
                     right_on='scaffold_id', 
                     how='left')
    
    # Reorder columns to have mag_id first
    bgc_df = bgc_df[['mag_id', 'scaffold_name', 'bgc_file', 'bgc_type', 'bgc_length']]
    
    bgc_df.to_csv(args.output_bgc_tsv, sep='\t', index=False)
    print(f"BGC data written to {args.output_bgc_tsv}")

    # Write peptide data to TSV if found
    if all_peptide_data:
        peptide_df = pd.DataFrame(all_peptide_data)
        # Merge with mag_scaffold mapping
        peptide_df = pd.merge(peptide_df, mag_scaffold_df,
                            left_on='scaffold_name',
                            right_on='scaffold_id',
                            how='left')
        peptide_df = peptide_df[['mag_id', 'scaffold_name', 'bgc_file', 'leader_seq', 'core_seq']]
        peptide_df.to_csv(args.output_peptide_tsv, sep='\t', index=False)
        print(f"Lanthipeptides found and written to {args.output_peptide_tsv}.")

        # Write FASTA records
        SeqIO.write(all_fasta_records, args.output_peptide_fasta, "fasta")
        print(f"Core peptide sequences written to {args.output_peptide_fasta}.")
    else:
        print(f"No lanthipeptides found in the input GBK files. No peptide files generated.")

if __name__ == "__main__":
    main()