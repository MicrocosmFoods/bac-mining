import argparse
import csv
import json
from Bio import SeqIO

def read_fasta(fasta_file):
    """Read a FASTA file using BioPython and return a dictionary of sequences."""
    sequences = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}
    return sequences

def extract_peptide_sequences(
    data,
    protein_fasta_file,
    proteins_output_file,
    protein_peptides_output_file,
    predictions_output_file,
    nucleotide_fasta_file=None,
    nucleotides_output_file=None,
    nucleotide_peptides_output_file=None,
):
    protein_sequences = read_fasta(protein_fasta_file)
    nucleotide_sequences = read_fasta(nucleotide_fasta_file) if nucleotide_fasta_file else {}

    protein_records = []
    nucleotide_records = []
    protein_peptide_records = []
    nucleotide_peptide_records = []
    predictions = []

    for protein_key, protein_info in data["PREDICTIONS"].items():
        protein_id = protein_key.split()[0][1:]
        peptides = protein_info.get("peptides")
        protein_sequence = protein_sequences.get(protein_id)
        nucleotide_sequence = nucleotide_sequences.get(protein_id)

        if protein_sequence and peptides:
            protein_records.append(SeqRecord(Seq(protein_sequence), id=protein_id, description=""))
            for peptide in peptides:
                start, end, peptide_class = peptide["start"], peptide["end"], peptide["type"]
                protein_peptide_sequence = protein_sequence[start - 1 : end]
                peptide_id = f"{protein_id}_start{start}_end{end}"
                description = " ".join(
                    [
                        f"{key}:{value}"
                        for key, value in {
                            "start": start,
                            "end": end,
                            "type": "cleavage",
                            "class": peptide_class,
                            "prediction_tool": "deeppeptide",
                        }.items()
                    ]
                )
                protein_peptide_records.append(
                    SeqRecord(Seq(protein_peptide_sequence), id=peptide_id, description=description)
                )
                predictions.append(
                    [peptide_id, start, end, "cleavage", peptide_class, "deeppeptide"]
                )

                if nucleotide_sequence:
                    nucleotide_peptide_sequence = nucleotide_sequence[(start - 1) * 3 : end * 3]
                    nucleotide_peptide_records.append(
                        SeqRecord(
                            Seq(nucleotide_peptide_sequence),
                            id=peptide_id,
                            description=description,
                        ),
                    )

            if nucleotide_sequence and peptides:
                nucleotide_records.append(
                    SeqRecord(Seq(nucleotide_sequence), id=protein_id, description="")
                )

    SeqIO.write(protein_records, proteins_output_file, "fasta")
    SeqIO.write(protein_peptide_records, protein_peptides_output_file, "fasta")
    if nucleotide_fasta_file:
        SeqIO.write(nucleotide_records, nucleotides_output_file, "fasta")
        SeqIO.write(nucleotide_peptide_records, nucleotide_peptides_output_file, "fasta")

    with open(predictions_output_file, "w", newline="") as predictions_out:
        writer = csv.writer(predictions_out, delimiter="\t")
        writer.writerow(
            ["peptide_id", "start", "end", "peptide_type", "peptide_class", "prediction_tool"]
        )
        writer.writerows(predictions)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract peptide sequences from DeepPeptide JSON.")
    parser.add_argument(
        "--json_file", type=str, required=True, help="The JSON file output by DeepPeptide."
    )
    parser.add_argument(
        "--protein_fasta_file",
        type=str,
        required=True,
        help="The protein FASTA file input to DeepPeptide.",
    )
    parser.add_argument(
        "--proteins_output_file", type=str, required=True, help="Output file path for proteins."
    )
    parser.add_argument(
        "--protein_peptides_output_file",
        type=str,
        required=True,
        help="Output file path for peptide sequences in amino acid format.",
    )
    parser.add_argument(
        "--predictions_output_file",
        type=str,
        required=True,
        help="Output file path for predictions.",
    )
    parser.add_argument(
        "--nucleotide_fasta_file", type=str, help="Optional nucleotide FASTA file for genes."
    )
    parser.add_argument(
        "--nucleotides_output_file",
        type=str,
        help="Optional output file path for nucleotide sequences.",
    )
    parser.add_argument(
        "--nucleotide_peptides_output_file",
        type=str,
        help="Optional output file path for peptide sequences in nucleotide format",
    )

    args = parser.parse_args()

    main(args)