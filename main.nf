#! /usr/bin/env nextflow

// Description
// Mine bacterial MAGs from fermented foods for peptides and BGCs

nextflow.enable.dsl=2

params.threads=4
params.outdir=null

log.info """\

MINE FERMENTED FOOD BACTERIAL MAGS FOR BIOACTIVE PEPTIDES AND BGCS
=================================================================
input_mags                      : $params.input_mags
outdir                          : $params.outdir
threads                         : $params.threads
"""

// define channels and workflow steps
mag_files = Channel.fromPath("${params.input_mags}/*.fa")
    .map { file -> 
        def baseName = file.getBaseName()
        return [file, baseName]
    }

workflow {
    // get small ORF predictions with smorfinder
    smorfinder(mag_files)
    smorf_proteins = smorfinder.out.faa_file
    // deepsig predictions
    deepsig(smorf_proteins)
    // peptides.py sequence characterization
    characterize_peptides(smorf_proteins)
    // DIAMOND seq similarity to Peptipedia peptide sequences of interest
    // input seqs to make comparisons, diamond db, diamond searches
    // autopeptideml tool for existing models of interest to compare? or do that outside of this?

}

process smorfinder {
    tag "${genome_name}_smorfinder"
    publishDir "${params.outdir}/smorfinder", mode: 'copy'

    conda: "envs/smorfinder.yml"

    input:
    tuple path(fasta), val(genome_name)

    output:
    tuple val(genome_name), path("*.gff"), emit: gff_file
    tuple val(genome_name), path("*.faa"), emit: faa_file
    tuple val(genome_name), path("*.ffn"), emit: ffn_file
    tuple val(genome_name), path("*.tsv"), emit: tsv_file

    script:
    """
    smorf single ${fasta} -o ${genome_name}
    ln -s ${genome_name}/${genome_name}.gff
    ln -s ${genome_name}/${genome_name}.faa
    ln -s ${genome_name}/${genome_name}.ffn
    ln -s ${genome_name}/${genome_name}.tsv
    """

}

process deepsig{
    tag "${genome_name}_deepsig"
    publishDir "${params.outdir}/deepsig", mode: 'copy'

    conda "envs/deepsig.yml"

    input: 
    tuple val(genome_name), path(faa_file)

    output: 
    path("*.tsv"), emit: deepsig_tsv

    script: 
    """
    deepsig -f ${faa_file} -o ${genome_name}.tsv -k gramp
    """

}

process characterize_peptides {
    tag "${genome_name}_characterize_peptides"
    publishDir "${params.outdir}/peptide_characterization", mode: 'copy'

    conda "envs/peptides.yml"

    input:
    tuple val(genome_name), path(faa_file)

    output: 
    path("*.tsv"), emit: peptides_tsv

    script:
    """
    python ${baseDir}/bin/characterize_peptides.py ${faa_file} ${genome_name}_peptide_characteristics.tsv
    """
}
