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
    
    // get small ORF predictions with smorfinder and combine into a single FASTA
    smorfinder(genome_fastas)
    smorf_proteins = smorfinder.out.faa_file.collect()
    combine_smorf_proteins(smorf_proteins)
    all_smorf_proteins = combine_smorf_proteins.out.combined_smorf_proteins

    // predict ORFs with pyrodigal
    pyrodigal(genome_fastas)
    predicted_orfs_gbks = pyrodigal.out.predicted_orfs_gbk
    predicted_orfs_proteins = pyrodigal.out.predicted_orfs_faa

    // predict encrypted peptides and combine into a single FASTA
    predict_encrypted_peptides(predicted_orfs_proteins)
    all_encrypted_peptides_fastas = predict_encrypted_peptides.out.encrypted_peptides_results_fasta.collect()
    combine_encrypted_peptides(all_encrypted_peptides_fastas)
    all_encrypted_peptides = combine_encrypted_peptides.out.combined_encrypted_peptides

    // predict cleavage peptides with deeppeptide, extract sequences from json, and combine into a single FASTA
    predict_cleavage_peptides(predicted_orfs_proteins)
    cleavage_peptides_json = predict_cleavage_peptides.out.cleavage_peptides_json
    cleavage_input_ch = cleavage_peptides_json.join(predicted_orfs_proteins, by: 0)
    extract_cleavage_peptides_json(cleavage_input_ch)
    all_cleavage_peptides_fastas = extract_cleavage_peptides_json.out.cleavage_peptides_fasta.collect()

    // antismash predictions
    antismash_input_ch = predicted_orfs_gbks.combine(antismash_db_ch)
    antismash(antismash_input_ch)
    antismash_gbk_files = antismash.out.gbk_results
    all_antismash_gbk_files = antismash_gbk_files.map{ it[1] }.collect()

    // extract antismash gbks into summary tsv, ripp peptides
    extract_gbks(all_antismash_gbk_files, genome_stb_tsv)
    all_core_ripp_peptides = extract_gbks.out.bgc_peptides_fasta

    // split out each peptide prediction tool into separate fastas by genome
    split_peptide_fastas_by_genome(all_smorf_proteins, all_encrypted_peptides, all_cleavage_peptides, all_core_ripp_peptides, genome_stb_tsv)
    per_genome_peptides_ch = split_peptide_fastas_by_genome.out.split_peptide_fastas
        .flatten()
        .map { file ->
            def genome_name = file.baseName
            return tuple(genome_name, file)
        }

    // cluster peptides per genome across prediction tools at 100% to get rid of redundancies
    mmseqs_100_cluster(per_genome_peptides_ch)
    all_nonredundant_peptides_ch = mmseqs_100_cluster.out.rep_seqs.collect()

    // summarize peptide and BGC counts per genome
    process_molecule_counts(all_nonredundant_peptides_ch, all_core_ripp_peptides, genome_stb_tsv)
    peptide_bgc_counts_summary = process_molecule_counts.out.peptide_bgc_counts_summary

    // run kofamscan annotations on all predicted proteins
    kofam_scan_annotation(predicted_orfs_proteins, kofam_db_ch)
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
    publishDir "${params.outdir}/main_results", mode: 'copy'

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
    path("*.csv"), emit: encrypted_peptides_results_csv
    path("*.fasta"), emit: encrypted_peptides_results_fasta

    script:
    """
    python ${baseDir}/bin/encrypted_peptides.py ${predicted_orfs_proteins} -o ${genome_name}_encrypted_peptides_results.csv -f ${genome_name}_encrypted_peptides_results.fasta
    """
}

process combine_encrypted_peptides {
    tag "combine_encrypted_peptides"
    publishDir "${params.outdir}/main_results", mode: 'copy'

    memory = "10 GB"
    cpus = 1

    container "quay.io/biocontainers/mulled-v2-949aaaddebd054dc6bded102520daff6f0f93ce6:aa2a3707bfa0550fee316844baba7752eaab7802-0"

    input:
    path(encrypted_peptides)

    output:
    path("*.fasta"), emit: combined_encrypted_peptides

    script:
    """
    python ${baseDir}/bin/combine_fastas.py ${encrypted_peptides.join(' ')} combined_encrypted_peptides.fasta
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
    tuple val(genome_name), path("*.json"), emit: cleavage_peptides_json

    script:
    """
    python3 predict.py --fastafile ${predicted_orfs_proteins} --output_dir ${genome_name} --output_fmt json
    mv ${genome_name}/*.json ./
    """
}

process extract_cleavage_peptides_json {
    tag "${genome_name}_extract_cleavage_peptides_json"
    publishDir "${params.outdir}/cleavage_peptides", mode: 'copy'

    memory = "10 GB"
    cpus = 1

    container "public.ecr.aws/biocontainers/mulled-v2-949aaaddebd054dc6bded102520daff6f0f93ce6:aa2a3707bfa0550fee316844baba7752eaab7802-0"

    input:
    tuple val(genome__name), path(deeppeptide_json), path(protein_faa)

    output:
    tuple val(genome_name), path("*_parent_proteins.faa"), emit: parent_proteins_faa
    tuple val(genome_name), path("*_peptides.faa"), emit: cleavage_peptides_fasta
    tuple val(genome_name), path("*.tsv"), emit: cleavage_peptides_tsv

    script:
    """
        python ${baseDir}/bin/extract_cleavage_peptides_json.py \
            --json_file ${deeppeptide_json} \
            --protein_fasta_file ${protein_faa} \
            --proteins_output_file ${genome_name}_parent_proteins.faa \
            --protein_peptides_output_file ${genome_name}_peptides.faa \
            --predictions_output_file ${genome_name}.tsv
    """
}

process combine_cleavage_peptides {
    tag "combine_cleavage_peptides"
    publishDir "${params.outdir}/main_results", mode: 'copy'

    memory = "10 GB"
    cpus = 1

    container "quay.io/biocontainers/mulled-v2-949aaaddebd054dc6bded102520daff6f0f93ce6:aa2a3707bfa0550fee316844baba7752eaab7802-0"

    input:
    path(cleavage_peptides)

    output:
    path("*.fasta"), emit: combined_cleavage_peptides

    script:
    """
    python ${baseDir}/bin/combine_fastas.py ${cleavage_peptides.join(' ')} combined_cleavage_peptides.fasta
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
    path(gbk_files)
    path(genome_stb_tsv)

    output:
    path("antismash_summary.tsv"), emit: bgc_summary_tsv
    path("antismash_peptides.tsv"), emit: bgc_peptides_tsv, optional: true
    path("antismash_peptides.fasta"), emit: bgc_peptides_fasta, optional: true

    script:
    """
    python ${baseDir}/bin/extract_antismash_gbks.py ${gbk_files.join(' ')} ${genome_stb_tsv} antismash_summary.tsv antismash_peptides.tsv antismash_peptides.fasta
    """

}

process split_peptide_fastas_by_genome {
    tag "split_peptide_fastas_by_genome"
    publishDir "${params.outdir}/split_peptide_fastas_by_genome", mode: 'copy'

    memory = "10 GB"
    cpus = 1

    container "quay.io/biocontainers/mulled-v2-949aaaddebd054dc6bded102520daff6f0f93ce6:aa2a3707bfa0550fee316844baba7752eaab7802-0"

    input:
    path(smorf_fasta)
    path(encrypted_fasta)
    path(cleavage_fasta)
    path(ripp_fasta)
    path(genome_stb)

    output:
    path("*.fasta"), emit: split_peptide_fastas

    script:
    """
    python ${baseDir}/bin/split_peptide_fastas_by_genome.py \\
    ${smorf_fasta} \\
    ${encrypted_fasta} \\
    ${cleavage_fasta} \\
    ${ripp_fasta} \\
    ${genome_stb} \\
    --outdir ./
    """
}

process mmseqs_100_cluster {
    tag "mmseqs_100_cluster"
    publishDir "${params.outdir}/mmseqs_100_cluster", mode: 'copy'

    memory = "10 GB"
    cpus = 1

    container "public.ecr.aws/biocontainers/mmseqs2:15.0--h5168794_0"

    input:
    tuple val(genome_name), path(split_peptide_fasta)

    output:
    path("*_rep_seq.fasta"), emit: rep_seqs
    path("*.tsv"), emit: clusters_tsv

    script:
    """
    mmseqs easy-cluster ${split_peptide_fasta} ${genome_name} --min-seq-id 1 -c 0.8 --threads ${task.cpus}
    """
}

process summarize_peptide_counts {
    tag "summarize_peptide_counts"
    publishDir "${params.outdir}/main_results", mode: 'copy'

    memory = "10 GB"
    cpus = 1

    container "quay.io/biocontainers/mulled-v2-949aaaddebd054dc6bded102520daff6f0f93ce6:aa2a3707bfa0550fee316844baba7752eaab7802-0"

    input:
    path("peptides/*")
    path(bgc_summary_tsv)
    path(genome_stb_tsv)

    output:
    path("*.tsv"), emit: peptide_counts_tsv

    script:
    """
    mkdir -p peptides
    python ${baseDir}/bin/process_molecule_counts.py \\
    --peptide-dir ./peptides \\
    --bgc-summary ${bgc_summary_tsv} \\
    --genome-stb ${genome_stb_tsv} \\
    --output peptide_bgc_counts_summary.tsv
    """
}

process kofamscan_annotation {
    tag "${genome_name}_kofam_scan_annotation"
    publishDir "${params.outdir}/kofam_scan_annotation", mode: 'copy'

    memory = "15 GB"
    cpus = 6

    container "public.ecr.aws/biocontainers/kofamscan:1.0.0--0"

    input:
    tuple val(genome_name), path(faa_file)
    path(kegg_db_dir)

    output:
    path("*.tsv"), emit: kofamscan_tsv

    script:
    """
    exec_annotation --format detail-tsv --ko-list ${kegg_db_dir}/ko_list --profile ${kegg_db_dir}/profiles --cpu ${task.cpus} -o peptides_kofamscan_annotation.tsv ${peptides_fasta}
    """
}