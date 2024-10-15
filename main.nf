#! /usr/bin/env nextflow

// Description
// Mine bacterial MAGs from fermented foods for peptides and BGCs

nextflow.enable.dsl=2

params.threads=4
params.outdir=null

log.info """\

MINE FERMENTED FOOD BACTERIAL MAGS FOR BIOACTIVE PEPTIDES AND BGCs.

NOTE: YOU MUST PRE-DOWNLOAD THE ANTISMASH DB AND PROVIDE THE PATH
WITH --antismash_db. THIS WORKFLOW DOES NOT SUPPORT DOWNLOADING THE DB
AUTOMATICALLY.
=================================================================
input_mags                      : $params.input_mags
antismash_db                    : $params.antismash_db
peptides_fasta                  : $params.peptides_fasta
outdir                          : $params.outdir
threads                         : $params.threads
"""

// define channels and workflow steps
mag_files = Channel.fromPath("${params.input_mags}/*.fa")
    .map { file -> 
        def baseName = file.getBaseName()
        return [file, baseName]
    }

antismash_db_ch = channel.fromPath(params.antismash_db, checkIfExists: true)
peptides_db_ch = channel.fromPath(params.peptides_fasta, checkIfExists: true)

workflow {
    // get small ORF predictions with smorfinder
    smorfinder(mag_files)
    smorf_proteins = smorfinder.out.faa_file

    // predict ORFs with pyrdogial
    pyrodigal(mag_files)
    predicted_proteins = pyrdogial.out.faa

    // antismash predictions and extract info from GBKs
    antismash(predicted_proteins, antismash_db_ch)
    antismash_gbk_files = antismash.out.gbk_results
    extract_antismash_info(antismash_gbk_files)

    // combine RIPP FASTAs with smORF FASTAs for each genome
    // get stats on peptide FASTAs, count for each genome, redundancy metrics

    // deepsig predictions
    deepsig(smorf_proteins)

    // peptides.py sequence characterization
    characterize_peptides(smorf_proteins)

    // DIAMOND seq similarity to Peptipedia peptide sequences of interest
    make_diamond_db(peptides_db_ch)
    peptides_dmnd_db = make_diamond_db.out.peptides_diamond_db
    diamond_blastp(smorf_proteins, peptides_dmnd_db)

}

process smorfinder {
    tag "${genome_name}_smorfinder"
    publishDir "${params.outdir}/smorfinder", mode: 'copy'

    conda "envs/smorfinder.yml"

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

process pyrodigal {
    tag "${genome_name}_pyrodigal"
    
    conda "envs/pyrodigal.yml"

    input:
    tuple path(fasta), val(genome_name)

    output:
    tuple val(genome_name), path("*.fna"), emit: fna
    tuple val(genome_name), path("*.faa"), emit: faa

    script:
    """
    pyrodigal \\
        -i ${fasta} \\
        -f "gbk" \\
        -o "${prefix}.${output_format}" \\
        -d ${genome_name}.fna \\
        -a ${genome_name}.faa
    """
}

process antismash {
    tag "${genome_name}_antismash"
    publishDir "${params.outdir}/antismash", mode: 'copy'

    conda "envs/antismashlite.yml"

    input:
    tuple val(genome_name), path(faa_file) // faa_file in gbk format
    path(databases)

    output: 
    tuple val(genome_name), path("${genome_name}/*.json") , emit: json_results
    tuple val(genome_name), path("${genome_name}/*.log") , emit: log
    tuple val(genome_name), path("${genome_name}/*region*.gbk") , optional: true, emit: gbk_results

    script: 
    """
    antismash \\
        -c $task.cpus \\
        --output-dir ${genome_name} \\
        --output-basename ${genome_name} \\
        --logfile ${genome_name}/${genome_name}.log \\
        --databases $databases \\
        ${faa_file}
    """
}

process extract_antismash_info {
    tag "${genome_name}_extract_antismash_info"
    publishDir "${params.outdir}/antismash_info", mode: 'copy'

    conda "envs/biopython.yml"

    input:
    tuple val(genome_name), path(gbk_files)

    output:
    tuple val(genome_name), path("*_antismash_summary.tsv"), emit: antismash_summary_tsv
    tuple val(genome_name), path("*_antismash_peptides.tsv"), emit: antismash_peptides_tsv
    tuple val(genome_name), path("*_antismash_peptides.fasta"), emit: antismash_peptides_fasta

    script:
    """
    python3 ${baseDir}/bin/extract_bgc_info_gbk.py ${genome_name} ${gbk_files.join(' ')} ${genome_name}_antismash_summary.tsv ${genome_name}_antismash_peptides.tsv ${genome_name}_antismash_peptides.fasta
    """
}

process deepsig {
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

process make_diamond_db {
    tag "make_diamond_db"

    conda "envs/diamond.yml"

    input:
    path(peptides_fasta)

    output:
    path("*.dmnd"), emit: peptides_diamond_db

    script:
    """
    diamond makedb --in ${peptides_fasta} -d peptides_db.dmnd
    """
}

process diamond_blastp {
    tag "diamond_blastp"
    publishDir "${params.outdir}/diamond_blastp", mode: 'copy'

    conda "envs/diamond.yml"

    input:
    path(combined_predicted_peptides_fasta)
    path(peptides_diamond_db)

    output:
    path("*.tsv"), emit: blastp_hits_tsv

    script:
    """
    diamond blastp -d ${peptides_diamond_db} -q ${combined_predicted_peptides_fasta} --header simple \\
    --outfmt 6 qseqid sseqid full_sseq pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore
    """
}
