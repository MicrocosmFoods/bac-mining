#!/usr/bin/env python

import argparse
import os
import pandas as pd
from Bio import SeqIO

class AntismashGBKParser:
    """Parse GBK file from antiSMASH to extract relevant BGC information."""

    def __init__(self, gbk_file):
        self.gbk_file = gbk_file
        self.records = list(SeqIO.parse(self.gbk_file, 'genbank'))

    def parse_records(self):
        """Parse all records in the GBK file."""
        bgc_data_list = []

        for record in self.records:
            scaffold_name = record.id
            bgc_type = self.get_bgc_type(record)
            bgc_length = len(record.seq)
            bgc_completeness = self.get_bgc_completeness(record)

            bgc_data_list.append({
                'scaffold_name': scaffold_name,
                'bgc_file': os.path.basename(self.gbk_file),
                'bgc_type': bgc_type,
                'bgc_length': bgc_length,
                'bgc_completeness': bgc_completeness,
            })

        return bgc_data_list

    def get_bgc_type(self, record):
        """Return the type of BGC region."""
        bgc_type_list = []
        for feature in record.features:
            if feature.type == 'region':
                bgc_type_list += feature.qualifiers.get('product', [])
        return '+'.join(set(bgc_type_list))

    def get_bgc_completeness(self, record):
        """Return the completeness of the BGC."""
        for feature in record.features:
            if feature.type == 'region':
                contig_edge = feature.qualifiers.get("contig_edge", ["complete"])[0]
                if contig_edge != "complete":
                    return "incomplete"
        return "complete"

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Parse antiSMASH GBK files to summarize BGCs.")
    parser.add_argument("input_gbk", nargs='+', help="Paths to the input GBK files")
    parser.add_argument("output_bgc_tsv", help="Path to the output TSV file for BGC information")
    return parser.parse_args()

def main():
    args = parse_args()

    all_bgc_data = []

    for gbk_file in args.input_gbk:
        parser = AntismashGBKParser(gbk_file)
        bgc_data = parser.parse_records()
        all_bgc_data.extend(bgc_data)

    # Write BGC data to TSV
    bgc_df = pd.DataFrame(all_bgc_data)
    bgc_df.to_csv(args.output_bgc_tsv, sep='\t', index=False)
    print(f"BGC data written to {args.output_bgc_tsv}")

if __name__ == "__main__":
    main()