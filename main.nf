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
    // get assembly stats with quast
    mag_stats = quast_stats(mag_files.map{ it[0] }.collect())
    
    // get ORFs, proteins, & GBK output with prodigal
    prodigal_predictions(mag_files)
    // get small ORF predictions with smorfinder
    smorfinder(mag_files)

}

process quast_stats {
    tag "quast_stats"
    publishDir "${params.outdir}/quast", mode: 'copy', pattern:"*.tsv"

    conda "envs/quast.yml"

    input: 
    path("*")

    output:
    path("*.tsv"), emit: quast_tsv

    script:
    """
    quast *.fa --output-dir QUAST -t 1
    ln -s QUAST/report.tsv
    ln -s QUAST/transposed_report.tsv
    """
}

process prodigal_predictions {
    tag "${genome_name}_prodigal"
    publishDir "${params.outdir}/prodigal", mode: 'copy'

    conda "envs/prodigal.yml"

    input:
    tuple path(fasta), val(genome_name)

    output:
    tuple val(genome_name), path("*.gbk"), emit: gbk_file
    tuple val(genome_name), path("*.faa"), emit: faa_file
    tuple val(genome_name), path("*.fna"), emit: fna_file

    script:
    """
    prodigal -i ${fasta} -o ${genome_name}.gbk -a ${genome_name}.faa -d ${genome_name}.fna
    """
}

process smorfinder {
    tag "${genome_name}_smorfinder"
    publishDir "${params.outdir}/smorfinder", mode: 'copy'

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