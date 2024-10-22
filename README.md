# Mining bacterial genomes for bioactive molecules

This repository contains the `bacMAGmining` Nextflow workflow for mining bacterial genomes (either metagenome-assembled genomes or isolate genomes) for bioactive peptides and biosynthetic gene clusters (BGCs). 

The main input is a directory of bacterial genomes in FASTA format. The workflow predicts small ORFs (smORFs) and BGCs, and outputs summary statistics of these predictions. Note that this workflow is highly configured for custom purposes, such as inputting a specific genome metadata TSV for joining metadata with main results files. You can see an example of the required metadata TSV file in `test_data/metadata/`. 

## Workflow Usage

The pipeline can be run with either conda or docker using the `-profile` flag. All input genomes in the input directory should end in `.fa`. 

Importantly, the workflow does not handle automatic downloading of databases, so these need to be prepared beforehand and input as parameters to the workflow. This includes the antiSMASH database that should be downloaded [according to the antiSMASH documentation](https://docs.antismash.secondarymetabolites.org/install/). Additionally the pipeline compares the predicted peptides to an input database of peptides of your choice. We run the pipeline with peptides from the [peptipedia database](https://app.peptipedia.cl/).

```
nextflow run main.nf \\
--input_genomes <INPUT_DIRECTORY> \\
--outdir <OUTPUT_DIRECTORY> \\
--genome_metadata <GENOME_METADATA_TSV> \\
--antismash_db <ANTISMASH_DB_DIR> \\
--peptides_fasta <PEPTIDES_FILE_FOR_COMPARISON> \\
--threads <THREADS> \\
-profile <docker|conda>
```