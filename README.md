# Mining bacterial genomes for bioactive molecules

This repository contains the `bacMAGmining` Nextflow workflow for mining bacterial genomes (either metagenome-assembled genomes or isolate genomes) for bioactive peptides and biosynthetic gene clusters (BGCs) and optionally performing functional annotation.

The main input is a directory of bacterial genomes in FASTA format. The workflow predicts small ORFs, cleavage peptides, and BGCs, which includes predicting RiPPs. Additionally the workflow provides optional functional annotation of whole proteomes using Kofamscan. The main output of the workflow are sets of FASTA files for each predicted peptide type, summaries of BGC types, and overall summaries of counts of each molecule type.  

## Workflow Usage

This pipeline can only be run with docker due to dependencies, and this is designated with the `-profile` flag. All input genomes in the input directory should end in `.fa`. 

Importantly, the workflow does not handle automatic downloading of databases, so these need to be prepared beforehand and input as parameters to the workflow. This includes the antiSMASH database that should be downloaded [according to the antiSMASH documentation](https://docs.antismash.secondarymetabolites.org/install/) and the Kofamscan database that should be downloaded from [here](https://www.genome.jp/kegg/rest/).

```
nextflow run main.nf \\
--input_genomes <INPUT_DIRECTORY> \\
--outdir <OUTPUT_DIRECTORY> \\
--antismash_db <ANTISMASH_DB_DIR> \\
--kofam_db <KO_DB_DIR> \\
--functional_annotation <true|false> \\
--threads <THREADS> \\
-profile <docker|conda>
```