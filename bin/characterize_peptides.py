import argparse
import csv

import peptides
from Bio import SeqIO

"""
This script is from the Arcadia-Science/peptigate pipeline for using the pepitdes.py package.
See https://github.com/Arcadia-Science/peptigate/blob/main/scripts/characterize_peptides.py for original source code.
"""


def characterize_peptides(input_file, output_file):
    """
    Processes a multi-FASTA file of peptides, calculates physicochemical properties and descriptors
    for each peptide, and writes the results to an output TSV file.

    This function uses the `peptides` library to calculate peptide properties like aliphatic index,
    boman index, charge, hydrophobicity, instability index, isoelectric point, molecular weight, and
    z-scales (lipophilicity, steric properties, electronic properties, etc.).
    It assumes default arguments for all peptide measurements as defined in the `peptides` library.
    For a comprehensive list of available measurements and their optional arguments, refer to the
    `peptides` library documentation: https://peptides.readthedocs.io.

    Parameters:
    - input_file (str): Path to the input FASTA file containing amino acid sequences of peptides.
    - output_file (str): Path to the output TSV file where the peptide properties will be written.

    Each row in the output TSV file includes the peptide ID, sequence, and calculated properties.

    Note: This function writes directly to the output file and does not return any value.
    """
    with open(output_file, "w", newline="") as out_file:
        tsv_writer = csv.writer(out_file, delimiter="\t")
        tsv_writer.writerow(
            [
                "peptide_id",
                "aliphatic_index",
                "boman_index",
                "charge",
                "hydrophobicity",
                "instability_index",
                "isoelectric_point",
                "molecular_weight",
                "pd1_residue_volume",
                "pd2_hydrophilicity",
                "z1_lipophilicity",
                "z2_steric_bulk_or_polarizability",
                "z3_polarity_or_charge",
                "z4_electronegativity_etc",
                "z5_electronegativity_etc",
            ]
        )

        for record in SeqIO.parse(input_file, "fasta"):
            peptide_sequence = peptides.Peptide(str(record.seq))
            aliphatic_index = peptide_sequence.aliphatic_index()
            boman_index = peptide_sequence.boman()
            charge = peptide_sequence.charge()
            hydrophobicity = peptide_sequence.hydrophobicity()
            instability_index = peptide_sequence.instability_index()
            isoelectric_point = peptide_sequence.isoelectric_point()
            molecular_weight = peptide_sequence.molecular_weight()
            physical_descriptors = peptide_sequence.physical_descriptors()
            zscales = peptide_sequence.z_scales()
            tsv_writer.writerow(
                [
                    record.id,
                    aliphatic_index,
                    boman_index,
                    charge,
                    hydrophobicity,
                    instability_index,
                    isoelectric_point,
                    molecular_weight,
                    physical_descriptors[0],
                    physical_descriptors[1],
                    zscales[0],
                    zscales[1],
                    zscales[2],
                    zscales[3],
                    zscales[4],
                ]
            )


def main():
    parser = argparse.ArgumentParser(description="Characterize peptides from a multi-fasta file.")
    parser.add_argument("input_file", type=str, help="Input multi-fasta file of amino acids")
    parser.add_argument("output_file", type=str, help="Output TSV file to write the results")

    args = parser.parse_args()

    characterize_peptides(args.input_file, args.output_file)


if __name__ == "__main__":
    main()