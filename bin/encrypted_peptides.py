###############################################################################
# Find Encrypted Peptides - This program searches a proteome to detect 
# encrypted antimicrobial peptides.
# Copyright (C) 2021 Marcelo Cardoso dos Reis Melo
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This script was modified from https://gitlab.com/machine-biology-group-public/amp_scoring/-/blob/main/FindEncryptedPeptides.py?ref_type=heads 
# by EAM for incorporation into the workflow.
###############################################################################

import argparse
from Bio import SeqIO
import numpy as np
import pandas as pd
import multiprocessing
import os
import time

def parse_arguments():
    parser = argparse.ArgumentParser(description='Search proteome for encrypted antimicrobial-like peptides')
    parser.add_argument('fasta_file', help='Path to input FASTA file')
    parser.add_argument('--output', '-o', default='encrypted_peptides_results.csv',
                      help='Output CSV file name (default: encrypted_peptides_results.csv)')
    parser.add_argument('--fasta-output', '-f', default='encrypted_peptides_results.fasta',
                        help='Output FASTA file name (default: encrypted_peptides_results.fasta)')
    parser.add_argument('--processes', '-p', type=int, default=20,
                      help='Number of CPU processes to use (default: 20)')
    parser.add_argument('--candidates', '-c', type=int, default=1000,
                      help='Number of candidates to save per peptide length (default: 1000)')
    parser.add_argument('--min-len', type=int, default=8,
                      help='Minimum peptide length (default: 8)')
    parser.add_argument('--max-len', type=int, default=50,
                      help='Maximum peptide length (default: 50)')
    return parser.parse_args()

# Hydriphobicity using the Parker-Gly0 scale
hScale = {
    "A": 0.229, "C": 0.274, "D": 0, "E": 0, "F": 0.949, "G": 0, "H": 0.229,
    "I": 0.873, "K": 0, "L": 0.949, "M": 0.631, "N": 0, "P": 0.229, "Q": 0,
    "R": 0.096, "S": 0, "T": 0.032, "V": 0.599, "W": 1, "Y": 0.484, "X": 0
}

def maxScorePaneJTB2017(L, m=0.9, n=1.1):
    """Calculates the maximum possible score for a peptide of a given length.
    This "ideal" peptide would have Z-many hydorphobic amino acids, and L-Z positively charged residues.
    Iterates of all possible values of "Z", the number of the most hydrophobic residue (with H = 1)
    Assumes all other residues are Arginines (positively charged, with "some" hydrophobicity).
    """
    return max([(nTer + cTer + L - Z)**m * ((L - Z)*hScale["R"] + Z)**n for Z in range(L)])

def calcPaneJTB2017(seq, m=0.9, n=1.1, maxScore=1):
    # Charges from N- and C-terminal groups
    charge = nTer + cTer
    # Count Arginines and Lysines
    charge += (seq.count("R") + seq.count("K")) * 1
    # Count Glutamic and Aspartic acids
    charge += (seq.count("E") + seq.count("D")) * -1
    
    # Hydrophobicity
    try:
        H = np.sum([hScale[aa] for aa in seq])
    except:
        H = 0
    
    if charge <= 0:
        score = 0
    else:
        score = (charge**m)*(H**n)/maxScore
    
    return score

def procPtn(inputDat):
    pep_set = set()
    ParM, ParN, pepLen, maxScore, recid, ptnSeq = inputDat
    ptnLen = len(ptnSeq)
    resList = []
    
    # Create list of peptides and scores
    for i in range(ptnLen - pepLen + 1):
        pepSeq = str(ptnSeq[i:i+pepLen])
        pep_set.add(pepSeq)
        score = calcPaneJTB2017(pepSeq, ParM, ParN, maxScore)
        resList.append([recid, i, pepLen, ptnLen, score, score*pepLen, pepSeq])
    
    # Avoid selecting overlapping candidate peptides
    startLimPrev = -1
    startLim = np.floor(pepLen/2)

    # Create a DF for all peptides in this protein
    tmpDFptn = pd.DataFrame.from_records(resList, columns=cols)
    startPos = []

    while startLimPrev < (ptnLen - pepLen - 1):
        tmpDFpi = tmpDFptn.loc[tmpDFptn.start <= startLim].loc[tmpDFptn.start > startLimPrev]
        newStart = tmpDFpi.loc[tmpDFpi['relativeScore'].idxmax()]["start"]
        startPos.append(newStart)

        startLimPrev = newStart + np.floor(pepLen/2)
        startLim = startLimPrev + np.floor(pepLen/2)
    
    # Select peptides based on starting position within protein sequence
    tmpDFptn = tmpDFptn.loc[tmpDFptn.start.isin(startPos)]
    # Exclude relative scores equal to zero
    tmpDFptn = tmpDFptn.loc[tmpDFptn.relativeScore > 0]
    # Remove repeated sequences (may happen in large proteins)
    tmpDFptn = tmpDFptn.drop_duplicates(subset="sequence", keep="last")
    
    return [tmpDFptn.copy(), pep_set]

def main():
    args = parse_arguments()

    genome_name = os.path.splitext(os.path.basename(args.fasta_file))[0]
    
    # Parameters for the function defined in Pane, JTB (2017)
    global ParM, ParN, nTer, cTer, cols
    ParM = 0.9
    ParN = 1.1
    nTer = 1.0
    cTer = 0.0
    cols = ["RecID", "start", "len", "ptnLen", "relativeScore", "absoluteScore", "sequence"]

    # Open and read FASTA file
    record_dict = SeqIO.to_dict(SeqIO.parse(args.fasta_file, "fasta"))
    print(f"There are {len(record_dict)} protein records in this database.")
    print("----")

    DFlist = []
    startTimeScan = time.time()

    # Initialize worker pool
    with multiprocessing.Pool(args.processes) as procPool:
        pep_set_all = 0
        startTimePepLen = time.time()
        
        # Loop over possible Encrypted AMP length
        for pepLen in range(args.min_len, args.max_len + 1):
            pepStartTime = time.time()
            print(f"Scanning for {pepLen}-residue long peptides.")
            
            maxScore = maxScorePaneJTB2017(pepLen, ParM, ParN)
            inputRecs = []
            
            for recid, recIter in record_dict.items():
                inputRecs.append([ParM, ParN, pepLen, maxScore, recid, str(recIter.seq)])
            
            ptnCounter = 0
            startTime = time.time()
            chunksize = 50
            
            # Use iterator-map to distribute records to workers
            imapIter = procPool.imap(procPtn, inputRecs, chunksize=chunksize)
            ptnsDFs = []
            
            for res, pep_set_ in imapIter:
                ptnsDFs.append(res)
                pep_set_all += len(pep_set_)
                ptnCounter += 1
                outputStride = (2*args.processes*chunksize)
                
                if not (ptnCounter % outputStride):
                    print(f"We have scanned {ptnCounter:>5d} peptides in {(time.time() - pepStartTime)/60:5.1f} minutes. "
                          f"The last {outputStride:>5d} peptides took {time.time() - startTime:5.1f} seconds", end="\r")
                    startTime = time.time()
            
            print("")  # Line break for next print
            
            # Combine results from all proteins
            ptnDFLen = pd.concat(ptnsDFs)
            
            # Remove repeated sequences (may happen in protein isoforms)
            ptnDFLen = ptnDFLen.drop_duplicates(subset="sequence", keep="first")
            
            # Sort and keep top candidates
            ptnDFLen = ptnDFLen.sort_values(by='relativeScore', ascending=False).iloc[0:args.candidates]

            print(f"{pepLen} - We found the highest relative score of {ptnDFLen['relativeScore'].max():.3f}. "
                  f"Elapsed Time: {time.time() - startTimePepLen:.2f}")
            
            startTimePepLen = time.time()
            DFlist.append(ptnDFLen.copy())
            
    # Combine all results and save
    allPepScanDF = pd.concat(DFlist)
    allPepScanDF.to_csv(args.output, index=False)

    # Create and output FASTA file of results
    with open(args.fasta_output, 'w') as fasta_output:
        for idx, row in allPepScanDF.iterrows():
            header = f">{genome_name}_encryptedpeptide_{idx+1} score={row['relativeScore']:.3f} length={row['len']}"
            sequence = row['sequence']
            fasta_output.write(f"{header}\n{sequence}\n")

    print(f"Total scan time was {time.time() - startTimeScan:.2f} seconds")
    print(f"Total number of peptides: {pep_set_all}")
    print(f"Results written to CSV: {args.output}")
    print(f"Results written to FASTA: {args.fasta_output}")

if __name__ == "__main__":
    main()