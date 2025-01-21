#! /usr/bin/env nextflow

// Description
// Mine bacterial MAGs from fermented foods for peptides and BGCs

nextflow.enable.dsl=2

params.threads=16
params.outdir=null

log.info """\

MINE FERMENTED FOOD BACTERIAL GENOMES FOR PEPTIDES AND BGCS, AND PERFORM
FUNCTIONAL ANNOTATION.

NOTE: YOU MUST PRE-DOWNLOAD THE ANTISMASH, PFAM, AND KOFAM DATABASES AND PROVIDE THE PATHS
WITH --antismash_db, --pfam_db, AND --kofam_db. THIS WORKFLOW DOES NOT SUPPORT DOWNLOADING 
DATABASES AUTOMATICALLY.
=================================================================
input_genomes                   : $params.input_genomes
genome_metadata                 : $params.genome_metadata
antismash_db                    : $params.antismash_db
pfam_db                         : $params.pfam_db
kofam_db                        : $params.kofam_db
outdir                          : $params.outdir
threads                         : $params.threads
"""

// define channels
// genome_fastas tuple with genome name and fasta filepath
genome_fastas = Channel.fromPath("${params.input_genomes}/*.fa")
    .map { file -> 
        def baseName = file.getBaseName()
        return [file, baseName]
    }

genome_metadata = channel.fromPath(params.genome_metadata)
antismash_db_ch = channel.fromPath(params.antismash_db)
kofam_db_ch = channel.fromPath(params.kofam_db)
pfam_db_ch = channel.fromPath(params.pfam_db)

// workflow steps
workflow {
    // make genome STB
    all_genome_fastas_ch = genome_fastas.map{ it[0] }.collect()
    make_genome_stb(all_genome_fastas_ch)
    genome_stb_tsv = make_genome_stb.out.stb_tsv
    
    // get small ORF predictions with smorfinder
    smorfinder(genome_fastas)
    smorf_proteins = smorfinder.out.faa_file

    // combine smorf proteins into a single FASTA
    combine_smorf_proteins(smorf_proteins.collect())
    combined_smorf_proteins = combine_smorf_proteins.out.combined_smorf_proteins

    // predict ORFs with pyrdogial
    pyrodigal(genome_fastas)
    predicted_orfs_gbks = pyrodigal.out.predicted_orfs_gbk
    predicted_orfs_proteins = pyrodigal.out.predicted_orfs_faa

    // predict encrypted peptides
    predict_encrypted_peptides(predicted_orfs_proteins)
    encrypted_peptides_results = predict_encrypted_peptides.out.encrypted_peptides_results

    // predict cleavage peptides with deeppeptide
    predict_cleavage_peptides(predicted_orfs_proteins)
    cleavage_peptides_outdir = predict_cleavage_peptides.out.cleavage_peptides_outdir

    // extract deeppeptide sequences from json
    extract_cleavage_peptides_json(cleavage_peptides_outdir)


    // cluster peptides at % identity
    mmseqs_95id_cluster(combined_smorf_proteins)
    nonredundant_smorfs = mmseqs_95id_cluster.out.nonredundant_seqs_fasta
    mmseqs_clusters = mmseqs_95id_cluster.out.cluster_summary_tsv

    // mmseqs cluster summaries and stats, merging with metadata
    summarize_mmseqs_clusters(mmseqs_clusters, nonredundant_smorfs, genome_metadata)

    // antismash predictions
    antismash_input_ch = predicted_orfs.combine(antismash_db_ch)
    antismash(antismash_input_ch)
    antismash_gbk_files = antismash.out.gbk_results
    all_antismash_gbk_files = antismash_gbk_files.map{ it[1] }.collect()

    // extract antismash gbks into summary tsv, ripp peptides
    extract_gbks(all_antismash_gbk_files, genome_stb_tsv)

    // run kofamscan annotations

}

process make_genome_stb {
    tag "make_genome_stb"
    publishDir "${params.outdir}/genomestb", mode: 'copy'

    memory = '10 GB'
    cpus = 1

    container "quay.io/biocontainers/mulled-v2-949aaaddebd054dc6bded102520daff6f0f93ce6:aa2a3707bfa0550fee316844baba7752eaab7802-0"

    input:
    path(fasta_files)

    output:
    path("*.tsv"), emit: stb_tsv

    script:
    """
    python ${baseDir}/bin/generate_genome_stb.py ${fasta_files.join(' ')} -o genomes_stb.tsv
    """

}

process smorfinder {
    tag "${genome_name}_smorfinder"
    publishDir "${params.outdir}/smorfinder", mode: 'copy'

    errorStrategy 'ignore'
    // rarely some genomes will fail for no discernible reason, skip over these

    memory = '10 GB'
    cpus = 1

    container "public.ecr.aws/v7p5x0i6/elizabethmcd/smorfinder:v0.2"

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

    memory = '10 GB'
    cpus = 1
    
    container "quay.io/biocontainers/mulled-v2-949aaaddebd054dc6bded102520daff6f0f93ce6:aa2a3707bfa0550fee316844baba7752eaab7802-0"

    input:
    path(smorf_proteins)

    output: 
    path("*.fasta"), emit: combined_smorf_proteins

    script:
    """
    python ${baseDir}/bin/combine_fastas.py ${smorf_proteins.join(' ')} combined_smorf_proteins.fasta
    """
    
}

process pyrodigal {
    tag "${genome_name}_pyrodigal"
    
    memory = "5 GB"
    cpus = 1

    container "public.ecr.aws/biocontainers/pyrodigal:3.4.1--py310h4b81fae_0"

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

process predict_encrypted_peptides {
    tag "${genome_name}_predict_encrypted_peptides"
    publishDir "${params.outdir}/encrypted_peptides", mode: 'copy'

    memory = "10 GB"
    cpus = 1

    container "public.ecr.aws/biocontainers/mulled-v2-949aaaddebd054dc6bded102520daff6f0f93ce6:aa2a3707bfa0550fee316844baba7752eaab7802-0"

    input:
    tuple val(genome_name), path(predicted_orfs_proteins)

    output:
    path("*.csv"), emit: encrypted_peptides_results

    script:
    """
    python ${baseDir}/bin/encrypted_peptides.py ${predicted_orfs_proteins} -o ${genome_name}_encrypted_peptides_results.csv
    """
}

process predict_cleavage_peptides {
    tag "${genome_name}_predict_cleavage_peptides"
    publishDir "${params.outdir}/cleavage_peptides", mode: 'copy'

    accelerator 1, type: 'nvidia-t4'
    cpus = 8

    container "public.ecr.aws/v7p5x0i6/elizabethmcd/deeppeptide:v0.1"

    input:
    tuple val(genome_name), path(predicted_orfs_proteins)

    output:
    path("*"), emit: cleavage_peptides_outdir

    script:
    """
    python3 predict.py --fastafile ${predicted_orfs_proteins} --output_dir ${genome_name} --output_fmt json
    """
}

process extract_cleavage_peptides_json {
    tag "${genome_name}_extract_cleavage_peptides_json"
    publishDir "${params.outdir}/cleavage_peptides", mode: 'copy'

    memory = "10 GB"
    cpus = 1


}

process mmseqs_95id_cluster {
    tag "mmseqs_95id_cluster"
    publishDir "${params.outdir}/mmseqs_95id_cluster", mode: 'copy'

    memory = '10 GB'
    cpus = 8
    
    container "public.ecr.aws/biocontainers/mmseqs2:15.6f452--pl5321h6a68c12_2"

    input:
    path(protein_fasta_file)
    
    output:
    path("*_rep_seq.fasta"), emit: nonredundant_seqs_fasta
    path("*.tsv"), emit: cluster_summary_tsv

    script:
    """
    mmseqs easy-cluster ${protein_fasta_file} nonredundant_smorf_proteins tmp --min-seq-id 0.95 --threads ${task.cpus}
    """   
}

process summarize_mmseqs_clusters {
    tag "summarize_mmseqs_clusters"
    publishDir "${params.outdir}/main_results/mmseqs_clusters", mode: 'copy'

    memory = "10 GB"
    cpus = 1

    container "quay.io/biocontainers/mulled-v2-949aaaddebd054dc6bded102520daff6f0f93ce6:aa2a3707bfa0550fee316844baba7752eaab7802-0"

    input:
    path(mmseqs_cluster_file)
    path(mmseqs_nonredundant_seqs)
    path(genome_metadata_tsv)

    output:
    path("mmseqs_summary.tsv"), emit: mmseqs_summary
    path("mmseqs_metadata.tsv"), emit: mmseqs_metadata
    path("mmseqs_substrate_counts.tsv"), emit: mmseqs_substrate_counts
    path("mmseqs_phylo_groups_counts.tsv"), emit: mmseqs_phylo_groups_counts

    script:
    """
    python ${baseDir}/bin/process_mmseqs_clusters.py ${mmseqs_cluster_file} ${mmseqs_nonredundant_seqs} ${genome_metadata_tsv} mmseqs_summary.tsv mmseqs_metadata.tsv mmseqs_substrate_counts.tsv mmseqs_phylo_groups_counts.tsv
    """ 
}

process antismash {
    tag "${genome_name}_antismash"
    publishDir "${params.outdir}/antismash", mode: 'copy'

    memory = "20 GB"
    cpus = 6

    container "public.ecr.aws/biocontainers/antismash-lite:7.1.0--pyhdfd78af_0"

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
        --genefinding-tool none \\
        ${gbk_file}
    """
}

process extract_gbks {
    tag "extract_gbks"
    publishDir "${params.outdir}/main_results/bgc_info", mode: 'copy'

    memory = "10 GB"
    cpus = 1

    container "quay.io/biocontainers/mulled-v2-949aaaddebd054dc6bded102520daff6f0f93ce6:aa2a3707bfa0550fee316844baba7752eaab7802-0"

    input: 
    path(gbk_files), path(mag_scaffold_tsv)

    output:
    path("antismash_summary.tsv"), emit: bgc_summary_tsv
    path("antismash_peptides.tsv"), emit: bgc_peptides_tsv, optional: true
    path("antismash_peptides.fasta"), emit: bgc_peptides_fasta, optional: true

    script:
    """
    python ${baseDir}/bin/extract_antismash_gbks.py ${gbk_files.join(' ')} ${mag_scaffold_tsv} antismash_summary.tsv antismash_peptides.tsv antismash_peptides.fasta
    """

}

process kofamscan_annotation {
    tag "${genome_name}_kofamscan_annotation"
    publishDir "${params.outdir}/kofamscan_annotation", mode: 'copy'

    memory = "10 GB"
    cpus = 1

    container "public.ecr.aws/biocontainers/kofamscan:1.3.0--pl5321h6a68c12_0"

    input:
    tuple val(genome_name), path(faa_file)

    output:
    tuple val(genome_name), path("*.tsv"), emit: kofamscan_annotation

    script:

}