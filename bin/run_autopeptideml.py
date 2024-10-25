import argparse
import os
from pathlib import Path

import pandas as pd
from autopeptideml.autopeptideml import AutoPeptideML
from autopeptideml.utils.embeddings import RepresentationEngine
from Bio import SeqIO

"""
This script is from the Arcadia-Science/peptigate pipeline for interacting with the AutoPeptideML package.
See https://github.com/Arcadia-Science/peptigate/blob/main/scripts/run_autopeptideml.py for original source code.
"""


def read_fasta(input_fasta):
    """
    Reads a FASTA file and returns a pandas DataFrame with IDs and sequences.

    Args:
    input_fasta (str): Path to the FASTA file.

    Returns:
    pd.DataFrame: DataFrame with columns 'ID' and 'sequence'.
    """
    sequences = []
    for seq_record in SeqIO.parse(input_fasta, "fasta"):
        sequences.append({"peptide_id": seq_record.id, "sequence": str(seq_record.seq)})
    return pd.DataFrame(sequences)


def predict_sequences(
    df, model_folder, model_name, threads=6, seed=42, batch_size=64, delete=True, tmp_dirname="tmp"
):
    """
    Predicts peptide sequence bioactivity using AutoPeptideML and returns the predictions DataFrame.

    Args:
    df (pd.DataFrame): DataFrame with sequences to predict.
    model_folder (str): Path to the model folder.
    model_name (str): Name of the model. Used to rename "prediction" column to output name.
    threads (int): Number of threads used to run the prediction.
    seed (int): Random seed.
    batch_size (int): Number of peptide sequences to compute in each batch.
    delete (log): Whether to delete the unmodified CSV file output by autopeptideml.
    tmp_dirname (str): Directory name supplied to AutoPeptideML's outputdir argument.

    Returns:
    pd.DataFrame: DataFrame with predictions.
    """
    autopeptideml = AutoPeptideML(verbose=True, threads=threads, seed=seed)
    representation_engine = RepresentationEngine(model="esm2-8m", batch_size=batch_size)

    predictions = autopeptideml.predict(
        df=df, re=representation_engine, ensemble_path=model_folder, outputdir=tmp_dirname
    )
    predictions.rename(columns={"prediction": model_name}, inplace=True)
    if delete:
        # autopeptideml writes prediction dataframe with uninformative column names to a
        # user-specified outputdir. This function specifies that folder as a temporary directory.
        # When delete == True, this function removes the output file written there.
        tmp_dirname = Path(tmp_dirname)
        os.remove(tmp_dirname / "predictions.csv")
    return predictions


def save_predictions(predictions, output_path):
    """
    Saves the predictions DataFrame to a TSV file.

    Args:
    predictions (pd.DataFrame): DataFrame with predictions produced by predicut_sequences().
    output_path (str): Path to save the TSV file.
    """
    predictions.to_csv(output_path, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser(description="Predict sequences using AutoPeptideML.")
    parser.add_argument("--input_fasta", required=True, help="Path to the FASTA file.")
    parser.add_argument("--model_folder", required=True, help="Path to the model folder.")
    parser.add_argument("--model_name", required=True, help="Name of the model.")
    parser.add_argument("--output_tsv", required=True, help="Path to the output TSV file.")

    args = parser.parse_args()

    df = read_fasta(args.input_fasta)
    predictions = predict_sequences(df, args.model_folder, args.model_name)
    save_predictions(predictions, args.output_tsv)


if __name__ == "__main__":
    main()