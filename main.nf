#! /usr/bin/env nextflow

// Description
// Mine bacterial MAGs from fermented foods for peptides and BGCs

nextflow.enable.dsl=2

params.threads=16
params.outdir=null

log.info """\

MINE FERMENTED FOOD BACTERIAL GENOMES FOR BIOACTIVE PEPTIDES AND BGCs.

NOTE: YOU MUST PRE-DOWNLOAD THE ANTISMASH DB AND PROVIDE THE PATH
WITH --antismash_db. THIS WORKFLOW DOES NOT SUPPORT DOWNLOADING THE DB
AUTOMATICALLY.
=================================================================
input_genomes                   : $params.input_genomes
antismash_db                    : $params.antismash_db
peptides_fasta                  : $params.peptides_fasta
outdir                          : $params.outdir
threads                         : $params.threads
"""

// define channels and workflow steps
genome_fastas = Channel.fromPath("${params.input_genomes}/*.fa")
    .map { file -> 
        def baseName = file.getBaseName()
        return [file, baseName]
    }

antismash_db_ch = channel.fromPath(params.antismash_db, checkIfExists: true)
peptides_db_ch = channel.fromPath(params.peptides_fasta, checkIfExists: true)

workflow {
    // get small ORF predictions with smorfinder
    smorfinder(genome_fastas)
    smorf_proteins = smorfinder.out.faa_file

    // combine smorf proteins into a single FASTA
    combine_smorf_proteins(smorf_proteins.collect())
    combined_smorf_proteins = combine_smorf_proteins.out.combined_smorf_proteins

    // cluster smorf proteins 100% identity and get representative seqs
    mmseqs_100id_cluster(combined_smorf_proteins)
    nonredundant_smorfs = mmseqs_100id_cluster.out.nonredundant_seqs_fasta

    // count total and nonredundant smORFs
    count_smorf_peptides(combined_smorf_proteins, nonredundant_smorfs)

    // predict ORFs with pyrdogial
    pyrodigal(genome_fastas)
    predicted_orfs = pyrodigal.out.predicted_orfs_gbk

    // antismash predictions and extract info from GBKs
    antismash_input_ch = predicted_orfs.combine(antismash_db_ch)
    antismash(antismash_input_ch)
    antismash_gbk_files = antismash.out.gbk_results
    extract_antismash_info(antismash_gbk_files)

    // deepsig predictions on combined, non-redundant smorf proteins
    // deepsig(nonredundant_smorfs)

    // peptides.py sequence characterization on combined, non-redundant smorf proteins
    characterize_peptides(nonredundant_smorfs)

    // DIAMOND seq similarity to Peptipedia peptide sequences of interest
    make_diamond_db(peptides_db_ch)
    peptides_dmnd_db = make_diamond_db.out.peptides_diamond_db
    diamond_blastp(nonredundant_smorfs, peptides_dmnd_db)

}

process smorfinder {
    tag "${genome_name}_smorfinder"
    publishDir "${params.outdir}/smorfinder", mode: 'copy'

    conda "envs/smorfinder.yml"

    input:
    tuple path(fasta), val(genome_name)

    output:
    path("*.gff"), emit: gff_file
    path("*.faa"), emit: faa_file
    path("*.ffn"), emit: ffn_file
    path("*.tsv"), emit: tsv_file

    script:
    """
    smorf single ${fasta} -o ${genome_name}
    ln -s ${genome_name}/${genome_name}.gff
    ln -s ${genome_name}/${genome_name}.faa
    ln -s ${genome_name}/${genome_name}.ffn
    ln -s ${genome_name}/${genome_name}.tsv
    """

}

process combine_smorf_proteins {
    tag "combine_smorf_proteins"
    publishDir "${params.outdir}/combined_smorf_proteins", mode: 'copy'

    conda "envs/biopython.yml"

    input:
    path(smorf_proteins)

    output: 
    path("*.fasta"), emit: combined_smorf_proteins

    script:
    """
    python ${baseDir}/bin/combine_fastas.py ${smorf_proteins.join(' ')} combined_smorf_proteins.fasta
    """
    
}

process mmseqs_100id_cluster {
    tag "mmseqs_100id_cluster"
    publishDir "${params.outdir}/mmseqs_100id_cluster", mode: 'copy'

    conda "envs/mmseqs2.yml"

    cpus = 8

    input:
    path(protein_fasta_file)
    
    output:
    path("*_rep_seq.fasta"), emit: nonredundant_seqs_fasta
    path("*.tsv"), emit: cluster_summary_tsv

    script:
    """
    mmseqs easy-cluster ${protein_fasta_file} nonredundant_smorf_proteins tmp --min-seq-id 1 --threads ${task.cpus}
    """   
}

process count_smorf_peptides {
    tag "count_smorf_peptides"
    publishDir "${params.outidr}/smorf_counts", mode: 'copy'

    conda "envs/biopython.yml"

    input:
    path(combined_smorf_proteins)
    path(nonredundant_smorf_proteins)

    output:
    path("*.tsv"), emit: smorf_counts_tsv

    script:
    
    """
    python ${baseDir}/bin/count_smorfs.py ${combined_smorf_proteins} ${nonredundant_smorf_proteins} genome_smorf_counts.tsv
    """
}

process pyrodigal {
    tag "${genome_name}_pyrodigal"
    
    conda "envs/pyrodigal.yml"

    input:
    tuple path(fasta), val(genome_name)

    output:
    tuple val(genome_name), path("*.fna"), emit: predicted_orfs_fna
    tuple val(genome_name), path("*.faa"), emit: predicted_orfs_faa
    tuple val(genome_name), path("*.gbk"), emit: predicted_orfs_gbk

    script:
    """
    pyrodigal \\
        -i ${fasta} \\
        -f "gbk" \\
        -o "${genome_name}.gbk" \\
        -d ${genome_name}.fna \\
        -a ${genome_name}.faa
    """
}

process antismash {
    tag "${genome_name}_antismash"
    publishDir "${params.outdir}/antismash", mode: 'copy'

    conda "envs/antismashlite.yml"

    cpus = 4

    input:
    tuple val(genome_name), path(gbk_file), path(databases)

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
        ${gbk_file}
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

    script:
    """
    python3 ${baseDir}/bin/extract_bgc_info_gbk.py ${gbk_files.join(' ')} ${genome_name}_antismash_summary.tsv
    """
}

process deepsig {
    tag "deepsig_predictions"
    publishDir "${params.outdir}/deepsig", mode: 'copy'

    conda "envs/deepsig.yml"

    accelerator 1, type: 'nvidia-t4'
    cpus = 8

    input: 
    path(faa_file)

    output: 
    path("*.tsv"), emit: deepsig_tsv

    script: 
    """
    deepsig -f ${faa_file} -o nonredundant_smorf_proteins_deepsig.tsv -k gramp -t ${task.cpus}
    """

}

process characterize_peptides {
    tag "characterize_peptides"
    publishDir "${params.outdir}/peptide_characterization", mode: 'copy'

    conda "envs/peptides.yml"

    input:
    path(faa_file)

    output: 
    path("*.tsv"), emit: peptides_tsv

    script:
    """
    python ${baseDir}/bin/characterize_peptides.py ${faa_file} nonredundant_smorf_proteins_peptide_characteristics.tsv
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
    path(faa_file)
    path(peptides_diamond_db)

    output:
    path("*.tsv"), emit: blastp_hits_tsv

    script:
    """
    diamond blastp -d ${peptides_diamond_db} -q ${faa_file} -o nonredundant_smorf_proteins_blast_results.tsv --header simple \\
    --outfmt 6 qseqid sseqid full_sseq pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore
    """
}
