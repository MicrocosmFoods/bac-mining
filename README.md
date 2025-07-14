# Mining bacterial genomes for bioactive molecules

This repository contains the `bacMAGmining` Nextflow workflow for mining bacterial genomes (either metagenome-assembled genomes or isolate genomes) for bioactive peptides and biosynthetic gene clusters (BGCs) and optionally performing functional annotation.

The main input is a directory of bacterial genomes in FASTA format. The workflow predicts small ORFs, cleavage peptides, and BGCs, which includes predicting RiPPs. Additionally the workflow provides optional functional annotation of whole proteomes using Kofamscan. The main output of the workflow are sets of FASTA files for each predicted peptide type, summaries of BGC types, and overall summaries of counts of each molecule type.  

## Workflow Usage

This pipeline can only be run with docker due to dependencies, and this is designated with the `-profile` flag. All input genomes in the input directory should end in `.fa`. 

Importantly, the workflow does not handle automatic downloading of databases, so these need to be prepared beforehand and input as parameters to the workflow. This includes the antiSMASH database that should be downloaded [according to the antiSMASH documentation](https://docs.antismash.secondarymetabolites.org/install/) and the Kofamscan database that should be downloaded from [here](https://www.genome.jp/kegg/rest/). You can optionally choose to run functional annotation with Kofamscan or not by setting the `functional_annotation` parameter to `true` or `false`. In either case you will still have to provide the path to the Kofamscan database.

```
nextflow run main.nf \\
--input_genomes <INPUT_DIRECTORY> \\
--outdir <OUTPUT_DIRECTORY> \\
--antismash_db <ANTISMASH_DB_DIR> \\
--kofam_db <KO_DB_DIR> \\
--functional_annotation <true|false> \\
--smorfinder_mode <single|pre_called> \\
--threads <THREADS> \\
-profile <docker|conda>
```

The Kofamscan database needs to be structured as follows:
```
profiles/
ko_list
```

Alternatively instead of downloading the entire Kofamscan database, you can download select KEGG HMM accessions (as long as you have checked they are in the ko_list file), and place them in the `profiles/` directory.

## Predicting smORF peptides with SmORFinder

We implement SmORFinder for predicting small orf (smORF) peptides. This program by default modifies a custom distribution of Prodigal that reduces the `--min-gene-length`. In our experience sometimes the custom binary of Prodigal within SmORFinder fails with known Prodgial bugs, however Prodigal is no longer supported and rather it is recommended to use Pyrodigal. Additionally, to layer on additional annotation types onto each genome such as KEGG Kofamscan, peptide information, BGCs, etc. the original genome needs to have ORFs called in a uniform way, and that isn't possible when running SmORFinder with their custom binary of Prodigal since the most updated Pyrodigal versions (even when running with the same flags) don't quite output in the same fashion.

To get around these problems, we have modified SmORFinder with an additional sub-workflow that can be called with the flag `smorfinder_mode pre_called` where you can pass pre-called ORFs using Pyrodigal and predict smORFs instead of using the custom Prodigal binary. In our experience this does sometimes lead to fewer numbers of smORFs predicted per genome, maybe 1 or 2 per genome since we don't implement the `--min-gene-length` or `--min-gene-overlap` flags in Pyrodigal. However this allows for different annotation types to be collated back together since the genome had ORFs predicted in a uniform fashion, and Pyrodigal is currently actively supported with bug fixes, whereas Prodigal is not.

## Check input genomes

Several steps in this workflow require scaffolds to be unique across all MAGs, or else those steps will fail. To check for duplicate scaffolds across MAGs and fix these by making them unique, a helper script is provided in the `scripts` directory. You can run it with: 

```
python scripts/check-mag-scaffolds.py \\
-i <INPUT_DIRECTORY> \\
-o <OUTPUT_DIRECTORY> \\
-d <DUPLICATES_FILE>
```

This will check for duplicate scaffolds and rename the offending duplicate scaffold by appending the filename of the genome (everything before .fa, such as the name of the MAG) to the scaffold ID. A list of the identified duplicates will be written to the file specified by `-d`, where the first column is the MAG, the second column is the original scaffold ID, and the third column is the new scaffold ID.

The workflow will also check for duplicate scaffold IDs when generating the genome STB file, and will raise an error if any are found. Importantly, the workflow will not modify the offending MAGs, but it will halt the workflow if any are found.

## Split input genomes into batches

In our experience, sometimes running the workflow on large amounts input genomes (~5,000-10,000+) can cause issues and the workflow will fail after long runtimes. To avoid this, you can split input genomes into batches and run the workflow on each batch separately. This can be done by providing a list of genomes names with the optional `genome_list` parameter which is a subset of genomes that are in the input directory. To make the batch lists, a helper script is provided in the `scripts` directory. You can run it with:

```
python scripts/prep-batch-genome-lists.py \\
-i <INPUT_DIRECTORY> \\
-b <BATCH_SIZE> \\
-o <OUTPUT_DIRECTORY>
```

Where `-b` is the maximum bathch size of genomes. In our experience running somewhere between 500-1000 genomes at a time should be fine. 

## Combine results from multiple runs

If you ran the workflow on multiple batches of genomes and want to combine the results into a single directory and set of results summary files, you can use the helper script `scripts/combine-batch-mining-results.py`. You can run it with: 

```
python scripts/combine-batch-mining-results.py \\
<INPUT_DIRECTORIES> \\
--output-dir <OUTPUT_DIRECTORY>
```

Where `<INPUT_DIRECTORIES>` is a list of directories containing the results from each batch run. The list of input directories should be separated by spaces, and needs to be the `main_results` folder downloaded from each batch run, where you have renamed each `main_results` folder from the name with some informative name and each ends in `-main-results` so you can track back to which batch those results came from. For example: 

```
python scripts/combine-batch-mining-results.py \\
2024-04-17_batch_1-main-results 2024-04-17_batch_2-main-results 2024-04-17_batch_3-main-results \\
--output-dir combined_results
```

This will combine the results into a single directory and set of results summary files. Each summary file will have a `batch_name` column that will have the name of the batch that the results came from, which is the name given to the `main_results` folder that is before the `-main-results` suffix.

## Create genome summaries

Currently the files in `main_results` are formatted where molecules of each type (such as smORFs, cleavage peptides, BGCs, KEGG Kofam annotations) are each separately summarized by all genomes in a single TSV file. And for all genomes and all molecule counts (excluding Kofam annotations) are summarized in a single file. 

If instead you want a TSV summary per genome of all annotations, you can do this with the helper script `scripts/create_genome_summaries.py`. 

To use this script you will need to point to the `main_results` directory and the directory where all the predicted ORFs in `.ffn` format are located (in the `predicted_orfs/`) subdirectory of the results output are all the predicted ORFs in .ffn, .faa, .gbk, and .gff formats. This script will output per genome a TSV summary of these annotations, where each row is a gene from the .ffn file and whether it has any of the annotations from the different programs. 