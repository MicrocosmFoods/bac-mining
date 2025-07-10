#!/usr/bin/env python

import argparse
import pandas as pd
import os
from Bio import SeqIO
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description="Combine all annotation information for a single genome")
    
    parser.add_argument("--ffn", required=True, help="Input .ffn file (nucleotide sequences)")
    parser.add_argument("--smorfinder-tsv", help="Smorfinder results TSV file")
    parser.add_argument("--deeppeptide-tsv", help="DeepPeptide results TSV file")
    parser.add_argument("--antismash-tsv", help="AntiSMASH results TSV file")
    parser.add_argument("--kofamscan-tsv", help="Kofamscan results TSV file")
    parser.add_argument("--output", required=True, help="Output TSV file path")
    
    return parser.parse_args()

def extract_locus_info_from_ffn(ffn_file):
    """Extract locus tag, start, and stop positions from .ffn file."""
    locus_info = {}
    
    for record in SeqIO.parse(ffn_file, "fasta"):
        # Parse the header to extract locus tag and position info
        # Example header: ">locus_tag_001 [location=1..150]"
        header = record.description
        
        # Extract locus tag (everything before the first space)
        locus_tag = record.id
        
        # Try to extract start and stop from location in description
        start_pos = "NA"
        stop_pos = "NA"
        
        if "[location=" in header:
            location_part = header.split("[location=")[1].split("]")[0]
            if ".." in location_part:
                start_pos, stop_pos = location_part.split("..")
                # Remove any strand indicators like "complement(" or ")"
                start_pos = start_pos.replace("complement(", "").replace("(", "")
                stop_pos = stop_pos.replace(")", "")
        
        locus_info[locus_tag] = {
            'start': start_pos,
            'stop': stop_pos
        }
    
    return locus_info

def load_smorfinder_results(tsv_file):
    """Load smorfinder results into a dictionary keyed by locus tag."""
    if not tsv_file or not os.path.exists(tsv_file):
        return {}
    
    results = {}
    try:
        df = pd.read_csv(tsv_file, sep='\t')
        for _, row in df.iterrows():
            locus_tag = row.get('locus_tag', '')
            if locus_tag:
                results[locus_tag] = {
                    'smorf_type': row.get('smorf_type', 'NA'),
                    'smorf_annotation': row.get('annotation', 'NA')
                }
    except Exception as e:
        print(f"Warning: Could not parse smorfinder file {tsv_file}: {e}")
    
    return results

def load_deeppeptide_results(tsv_file):
    """Load DeepPeptide results into a dictionary keyed by locus tag."""
    if not tsv_file or not os.path.exists(tsv_file):
        return {}
    
    results = {}
    try:
        df = pd.read_csv(tsv_file, sep='\t')
        for _, row in df.iterrows():
            locus_tag = row.get('locus_tag', '')
            if locus_tag:
                results[locus_tag] = {
                    'peptide_class': row.get('peptide_class', 'NA'),
                    'peptide_annotation': row.get('annotation', 'NA')
                }
    except Exception as e:
        print(f"Warning: Could not parse deeppeptide file {tsv_file}: {e}")
    
    return results

def load_antismash_results(tsv_file):
    """Load AntiSMASH results into a dictionary keyed by locus tag."""
    if not tsv_file or not os.path.exists(tsv_file):
        return {}
    
    results = {}
    try:
        df = pd.read_csv(tsv_file, sep='\t')
        for _, row in df.iterrows():
            locus_tag = row.get('gene_name', '')
            if locus_tag:
                results[locus_tag] = {
                    'bgc_type': row.get('bgc_type', 'NA'),
                    'bgc_annotation': row.get('gene_function', 'NA')
                }
    except Exception as e:
        print(f"Warning: Could not parse antismash file {tsv_file}: {e}")
    
    return results

def load_kofamscan_results(tsv_file):
    """Load Kofamscan results into a dictionary keyed by locus tag."""
    if not tsv_file or not os.path.exists(tsv_file):
        return {}
    
    results = {}
    try:
        df = pd.read_csv(tsv_file, sep='\t')
        for _, row in df.iterrows():
            locus_tag = row.get('gene_name', '')
            if locus_tag:
                results[locus_tag] = {
                    'kegg_ko': row.get('ko', 'NA'),
                    'kegg_annotation': row.get('definition', 'NA')
                }
    except Exception as e:
        print(f"Warning: Could not parse kofamscan file {tsv_file}: {e}")
    
    return results

def main():
    args = parse_args()
    
    # Extract locus information from .ffn file
    print(f"Extracting locus information from {args.ffn}")
    locus_info = extract_locus_info_from_ffn(args.ffn)
    
    # Load all annotation results
    print("Loading annotation results...")
    smorf_results = load_smorfinder_results(args.smorfinder_tsv)
    deep_results = load_deeppeptide_results(args.deeppeptide_tsv)
    antismash_results = load_antismash_results(args.antismash_tsv)
    kofam_results = load_kofamscan_results(args.kofamscan_tsv)
    
    # Create comprehensive summary
    print("Creating comprehensive summary...")
    summary_rows = []
    
    for locus_tag in sorted(locus_info.keys()):
        row = {
            'locus_tag': locus_tag,
            'start': locus_info[locus_tag]['start'],
            'stop': locus_info[locus_tag]['stop']
        }
        
        # Add KEGG annotation
        if locus_tag in kofam_results:
            row['kegg_ko'] = kofam_results[locus_tag]['kegg_ko']
            row['kegg_annotation'] = kofam_results[locus_tag]['kegg_annotation']
        else:
            row['kegg_ko'] = 'NA'
            row['kegg_annotation'] = 'NA'
        
        # Add AntiSMASH BGC information
        if locus_tag in antismash_results:
            row['bgc_type'] = antismash_results[locus_tag]['bgc_type']
            row['bgc_annotation'] = antismash_results[locus_tag]['bgc_annotation']
        else:
            row['bgc_type'] = 'NA'
            row['bgc_annotation'] = 'NA'
        
        # Add peptide annotations
        peptide_types = []
        peptide_annotations = []
        
        if locus_tag in smorf_results:
            peptide_types.append('smorf')
            peptide_annotations.append(smorf_results[locus_tag]['smorf_annotation'])
        
        if locus_tag in deep_results:
            peptide_types.append(deep_results[locus_tag]['peptide_class'])
            peptide_annotations.append(deep_results[locus_tag]['peptide_annotation'])
        
        if peptide_types:
            row['peptide_types'] = ';'.join(peptide_types)
            row['peptide_annotations'] = ';'.join(peptide_annotations)
        else:
            row['peptide_types'] = 'NA'
            row['peptide_annotations'] = 'NA'
        
        summary_rows.append(row)
    
    # Create and save summary dataframe
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(args.output, sep='\t', index=False)
    print(f"Summary saved to {args.output}")
    print(f"Total loci processed: {len(summary_rows)}")

if __name__ == "__main__":
    main() 